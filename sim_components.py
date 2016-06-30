import simpy
import random
import copy
from simpy.core import BoundClass
from simpy.resources import base
from heapq import heappush, heappop
from patch import *

class Packet(object):
  def __init__(self, time, size, _id, src="a", dst="z", flow_id=0):
    self.time = time
    self.size = size
    self._id = _id
    self.src = src
    self.dst = dst
    self.flow_id = flow_id
    self.ref_time = 0 # for casual use
    # for FJ and MDS Q implementation
    self.prev_hop_id = None
    self.entrance_time = 0
    self.job_id = None
  
  def deep_copy(self):
    p = Packet(time=self.time, size=self.size, _id=self._id, src=self.src, dst=self.dst, flow_id=self.flow_id)
    p.ref_time = self.ref_time
    p.prev_hop_id = self.prev_hop_id
    p.entrance_time = self.entrance_time
    p.job_id = self.job_id
    return p
  
  def __repr__(self):
    return "Packet[_id: {}, prev_hop_id: {}, job_id: {}]".\
      format(self._id, self.prev_hop_id, self.job_id)

class CPacket(object): # Control
  def __init__(self, _id):
    self._id = _id
  
  def __repr__(self):
    return "CPacket[_id= {}]".format(self._id)

class PacketGenerator(object):
  def __init__(self, env, _id, adist, sdist, initial_delay=0, finish=float("inf"), flow_id=0):
    self._id = _id
    self.env = env
    self.adist = adist
    self.sdist = sdist
    self.initial_delay = initial_delay
    self.finish = finish
    self.n_sent = 0
    self.flow_id = flow_id
    
    self.out = None
    self.action = env.process(self.run())  # starts the run() method as a SimPy process

  def run(self):
    yield self.env.timeout(self.initial_delay)
    while self.env.now < self.finish:
      # wait for next transmission
      yield self.env.timeout(self.adist() )
      self.n_sent += 1
      p = Packet(time=self.env.now, size=self.sdist(), _id=self.n_sent, src=self._id, flow_id=self.flow_id)
      self.out.put(p)

class S1_Q(object): # Memoryless service, 1 server
  def __init__(self, _id, env, serv_dist=None, rate=None, qlimit_n=None, qlimit_B=None, debug=False):
    self._id = _id
    self.env = env
    self.serv_dist = serv_dist
    self.rate = rate
    self.qlimit_n = qlimit_n
    self.qlimit_B = qlimit_B
    self.debug = debug
    
    self.p_list = []
    self.p_in_serv = None
    self.cancel = None
    self.n_recved = 0
    self.n_dropped = 0
    self.size_n = 0  # Current size of the queue in n
    self.size_B = 0  # Current size of the queue in bytes
    self.busy = 0  # Used to track if a packet is currently being sent
    self.wt_list = []
    self.qt_list = []
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.syncer = simpy.Store(env) # simpy.Resource(env, capacity=1)
    self.out = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    self.action = env.process(self.run_helper() )
  
  def __repr__(self):
    # return "S1_Q[_id= {}, rate= {}, qlimit_n= {}, qlimit_B= {}]".format(_id, self.rate, self.qlimit_n, self.qlimit_B)
    return "S1_Q[_id= {}]".format(self._id)
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      self.p_list.append(p)
      self.syncer.put(1)
  
  def run_helper(self): # To implement del from self.store
    while True:
      (yield self.syncer.get() )
      if len(self.p_list) == 0:
        continue # log(ERROR, "self.p_list is empty!") # May happen because of task cancellations
      self.p_in_serv = self.p_list.pop(0)
      self.wt_list.append(self.env.now - self.p_in_serv.ref_time)
      self.size_n -= 1
      self.size_B -= self.p_in_serv.size
      self.busy = 1
      if self.cancel is None:
        self.cancel = self.env.event()
      if self.serv_dist is None:
        if self.rate is None:
          log(ERROR, "self.serv_dist is None but self.rate is None too!")
          return 1
        yield (self.cancel | self.env.timeout(self.p_in_serv.size/self.rate) ) # service
      else:
        yield (self.cancel | self.env.timeout(self.serv_dist() ) ) # service
      if self.cancel is None: # task got cancelled
        sim_log(DEBUG, self.env, self, "cancelling", self.p_in_serv)
      else:
        self.qt_list.append(self.env.now - self.p_in_serv.ref_time)
        if self.out is not None:
          sim_log(DEBUG, self.env, self, "finished serv, forwarding", self.p_in_serv)
          self.p_in_serv.prev_hop_id = self._id
          self.out.put(self.p_in_serv)
      self.busy = 0
  
  def put(self, p):
    self.n_recved += 1
    p.ref_time = self.env.now
    sim_log(DEBUG, self.env, self, "recved", p)
    t_size_n = self.size_n + 1
    t_size_B = self.size_B + p.size
    if (self.qlimit_n is not None and t_size_n > self.qlimit_n) or \
       (self.qlimit_B is not None and t_size_B > self.qlimit_B):
      sim_log(DEBUG, self.env, self, "dropping", p)
      self.n_dropped += 1
      return
    else:
      self.size_n = t_size_n
      self.size_B = t_size_B
      return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      if self.p_in_serv.job_id == cp._id:
        if self.cancel is None:
          log(ERROR, "self.cancel is None!")
          return 1
        self.cancel.succeed()
        self.cancel = None
      
      for p in self.p_list:
        if p.job_id == cp._id:
          self.p_list.remove(p)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)
    
class QMonitor(object):
  def __init__(self, env, q, dist):
    self.q = q
    self.env = env
    self.dist = dist
    
    self.t_list = [] # time steps that the numbers polled from the q
    self.n_list = [] # num of packets in the q
    self.action = env.process(self.run() )
  
  def run(self):
    while True:
      yield self.env.timeout(self.dist() )
      
      self.t_list.append(self.env.now)
      self.n_list.append(len(self.q.store.items) + self.q.busy)

# *******************************************  MDS  ********************************************** #
class JSink(object): # Join
  def __init__(self, _id, env):
    self._id = _id
    self.env = env
    
    self.st_list = [] # total time in system
    
    self.store = simpy.Store(env)
    self.out = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
  
  def __repr__(self):
    return "JSink[_id= {}]".format(self._id)
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      self.st_list.append(self.env.now - p.entrance_time)
      if self.out is not None:
        self.out.put(p)
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    return self.store.put(p)

class JQ(object): # JoinQ for MDS; completion of any k tasks out of n means job comletion
  def __init__(self, env, k, input_qid_list):
    self.env = env
    self.k = k
    self.input_qid_list = input_qid_list
    
    self.input_id__pq_map = {i: [] for i in input_qid_list}
    self.qt_list = []
    self.length = 0 # maximum of the lengths of all pq's
    
    self.store = simpy.Store(env)
    self.out = None
    self.out_c = None
    self.action = env.process(self.run() )  # starts the run() method as a SimPy process
  
  def __repr__(self):
    return "JQ[k= {}, input_qid_list= [{}]]".format(self.k, ", ".join(self.input_qid_list) )
    
  def check_for_job_completion(self):
    now = self.env.now
    len_list = [len(pq) for i, pq in self.input_id__pq_map.items() ]
    num_non_zero = len([l for l in len_list if l > 0] )
    if num_non_zero > self.k:
      log(ERROR, "num_non_zero= {} > k= {}".format(num_non_zero, self.k) )
    elif num_non_zero < self.k:
      return 0
    
    ref_p = None
    for j, pq in self.input_id__pq_map.items():
      if len(pq) == 0:
        continue
      p = pq.pop(0)
      if ref_p is None:
        ref_p = p
      else:
        if (p.prev_hop_id == ref_p.prev_hop_id) or \
           (p.entrance_time != ref_p.entrance_time) or \
           (p.job_id != ref_p.job_id):
          log(ERROR, "supposed to be tasks of the same job;\np= {}\nref_p= {}".format(p, ref_p) )
          return 1
      self.qt_list.append(now - p.ref_time)
    # signal cancellation of the remaining tasks of the completed job
    self.out_c.put_c(CPacket(ref_p.job_id) )
    
    self.out.put(ref_p)
    self.length = max([len(pq) for i, pq in self.input_id__pq_map.items() ] )
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      if p.prev_hop_id not in self.input_qid_list:
        log(ERROR, "packet can NOT continue {}; packet= {}".format(self, p) )
        return 1
      self.input_id__pq_map[p.prev_hop_id].append(p)
      self.check_for_job_completion()
  
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.ref_time = self.env.now
    return self.store.put(p)

class MDSQ(object):
  def __init__(self, _id, env, k, qid_list, qserv_dist_list):
    self._id = _id
    self.env = env
    self.k = k
    self.qid_list = qid_list
    self.qserv_dist_list = qserv_dist_list
    
    self.num_q = len(qid_list)
    self.join_sink = JSink(_id, env)
    self.join_q = JQ(env, k, qid_list)
    self.join_q.out = self.join_sink
    self.join_q.out_c = self
    self.id_q_map = {}
    for i in range(self.num_q):
      qid = qid_list[i]
      m1_q = S1_Q(_id=qid, env=env, serv_dist=qserv_dist_list[i] )
      m1_q.out = self.join_q
      self.id_q_map[qid] = m1_q
    
    self.store = simpy.Store(env)
    self.store_c = simpy.Store(env)
    self.out = None
    self.action = env.process(self.run() ) # starts the run() method as a SimPy process
    self.action = env.process(self.run_c() )
    
    self.job_id_counter = 0
  
  def __repr__(self):
    return "MDSQ[k= {}, qid_list= [{}] ]".format(self.k, ", ".join(self.qid_list) )
  
  def run(self):
    while True:
      p = (yield self.store.get() )
      for i, q in self.id_q_map.items():
        q.put(p.deep_copy() )
      
  def put(self, p):
    sim_log(DEBUG, self.env, self, "recved", p)
    p.entrance_time = self.env.now
    self.job_id_counter += 1
    p.job_id = self.job_id_counter
    return self.store.put(p)
  
  def run_c(self):
    while True:
      cp = (yield self.store_c.get() )
      for i, q in self.id_q_map.items():
        q.put_c(cp)
  
  def put_c(self, cp):
    sim_log(DEBUG, self.env, self, "recved", cp)
    return self.store_c.put(cp)
  
# ################################################################################################ #
class RandomBrancher(object):
  """ A demultiplexing element that chooses the output port at random.

    Contains a list of output ports of the same length as the probability list
    in the constructor.  Use these to connect to other network elements.

    Parameters
    ----------
    env : simpy.Environment
      the simulation environment
    probs : List
      list of probabilities for the corresponding output ports
  """
  def __init__(self, env, probs):
    self.env = env

    self.probs = probs
    self.ranges = [sum(probs[0:n+1]) for n in range(len(probs))]  # Partial sums of probs
    if self.ranges[-1] - 1.0 > 1.0e-6:
      raise Exception("Probabilities must sum to 1.0")
    self.n_ports = len(self.probs)
    self.outs = [None for i in range(self.n_ports)]  # Create and initialize output ports
    self.n_recved = 0

  def put(self, pkt):
    self.n_recved += 1
    rand = random.random()
    for i in range(self.n_ports):
      if rand < self.ranges[i]:
        if self.outs[i]:  #  A check to make sure the output has been assigned before we put to it
          self.outs[i].put(pkt)
        return


class FlowDemux(object):
    """ A demultiplexing element that splits packet streams by flow_id.

    Contains a list of output ports of the same length as the probability list
    in the constructor.  Use these to connect to other network elements.

    Parameters
    ----------
    env : simpy.Environment
      the simulation environment
    outs : List
      list of probabilities for the corresponding output ports
  """
    def __init__(self, outs=None, default=None):
      self.outs = outs
      self.default = default
      self.n_recved = 0

    def put(self, pkt):
      self.n_recved += 1
      flow_id = pkt.flow_id
      if flow_id < len(self.outs):
        self.outs[flow_id].put(pkt)
      else:
        if self.default:
          self.default.put(pkt)

class SnoopSplitter(object):
  """ A snoop port like splitter. Sends the original packet out port 1
    and sends a copy of the packet out port 2.

    You need to set the values of out1 and out2.
  """
  def __init__(self):
    self.out1 = None
    self.out2 = None

  def put(self, pkt):
    pkt2 = copy.copy(pkt)
    if self.out1:
      self.out1.put(pkt)
    if self.out2:
      self.out2.put(pkt2)

"""
  Trying to implement a stamped/ordered version of the Simpy Store class.
  The "stamp" is used to sort the elements for removal ordering. This
  can be used in the implementation of sophisticated queueing disciplines, but
  would be overkill for fixed priority schemes.

  Copyright 2014 Greg M. Bernstein
  Released under the MIT license
"""

class StampedStorePut(base.Put):
  """ Put *item* into the store if possible or wait until it is.
    The item must be a tuple (stamp, contents) where the stamp is used to sort
    the content in the StampedStore.
  """
  def __init__(self, resource, item):
    self.item = item
    """The item to put into the store."""
    super(StampedStorePut, self).__init__(resource)


class StampedStoreGet(base.Get):
  """Get an item from the store or wait until one is available."""
  pass


class StampedStore(base.BaseResource):
  """Models the production and consumption of concrete Python objects.

  Items put into the store can be of any type.  By default, they are put and
  retrieved from the store in a first-in first-out order.

  The *env* parameter is the :class:`~simpy.core.Environment` instance the
  container is bound to.

  The *capacity* defines the size of the Store and must be a positive number
  (> 0). By default, a Store is of unlimited size. A :exc:`ValueError` is
  raised if the value is negative.

  """
  def __init__(self, env, capacity=float('inf')):
    super(StampedStore, self).__init__(env)
    if capacity <= 0:
      raise ValueError('"capacity" must be > 0.')
    self._capacity = capacity
    self.items = []  # we are keeping items sorted by stamp
    """List of the items within the store."""

  @property
  def capacity(self):
    """The maximum capacity of the store."""
    return self._capacity

  put = BoundClass(StampedStorePut)
  """Create a new :class:`StorePut` event."""

  get = BoundClass(StampedStoreGet)
  """Create a new :class:`StoreGet` event."""

  # We assume the item is a tuple: (stamp, packet). The stamp is used to
  # sort the packet in the heap.
  def _do_put(self, event):
    if len(self.items) < self._capacity:
      heappush(self.items, event.item)
      event.succeed()

  # When we return an item from the stamped store we do not
  # return the stamp but only the content portion.
  def _do_get(self, event):
    if self.items:
      event.succeed(heappop(self.items)[1])

"""
  A Set of components to enable simulation of various networking QoS scenarios.

  Copyright 2014 Dr. Greg Bernstein
  Released under the MIT license
"""


class ShaperTokenBucket(object):
  """ Models an ideal token bucket shaper. Note the token bucket size should be greater than the
    size of the largest packet that can occur on input. If this is not the case we always accumulate
    enough tokens to let the current packet pass based on the average rate. This maynot be
    the behavior you desire.

    Parameters
    ----------
    env : simpy.Environment
      the simulation environment
    rate : float
      the token arrival rate in bits
    b_size : Number
      a token bucket size in bytes
    peak : Number or None for infinite peak
      the peak sending rate of the buffer (quickest time two packets could be sent)

  """
  def __init__(self, env, rate, b_size, peak=None, debug=False):
    self.store = simpy.Store(env)
    self.rate = rate
    self.env = env
    self.out = None
    self.n_recved = 0
    self.n_sent = 0
    self.b_size = b_size
    self.peak = peak

    self.current_bucket = b_size  # Current size of the bucket in bytes
    self.update_time = 0.0  # Last time the bucket was updated
    self.debug = debug
    self.busy = 0  # Used to track if a packet is currently being sent ?
    self.action = env.process(self.run())  # starts the run() method as a SimPy process

  def run(self):
    while True:
      msg = (yield self.store.get())
      now = self.env.now
      #  Add tokens to bucket based on current time
      self.current_bucket = min(self.b_size, self.current_bucket + self.rate*(now-self.update_time)/8.0)
      self.update_time = now
      #  Check if there are enough tokens to allow packet to be sent
      #  If not we will wait to accumulate enough tokens to let this packet pass
      #  regardless of the bucket size.
      if msg.size > self.current_bucket:  # Need to wait for bucket to fill before sending
        yield self.env.timeout((msg.size - self.current_bucket)*8.0/self.rate)
        self.current_bucket = 0.0
        self.update_time = self.env.now
      else:
        self.current_bucket -= msg.size
        self.update_time = self.env.now
      # Send packet
      if not self.peak:  # Infinite peak rate
        self.out.put(msg)
      else:
        yield self.env.timeout(msg.size*8.0/self.peak)
        self.out.put(msg)
      self.n_sent += 1
      if self.debug:
        print(msg)

  def put(self, pkt):
    self.n_recved += 1
    return self.store.put(pkt)


class VirtualClockServer(object):
  """ Models a virtual clock server. For theory and implementation see:
    L. Zhang, Virtual clock: A new traffic control algorithm for packet switching networks,
    in ACM SIGCOMM Computer Communication Review, 1990, vol. 20, pp. 19.


    Parameters
    ----------
    env : simpy.Environment
      the simulation environment
    rate : float
      the bit rate of the port
    vticks : A list
      list of the vtick parameters (for each possible packet flow_id). We assume a simple assignment of
      flow id to vticks, i.e., flow_id = 0 corresponds to vticks[0], etc... We assume that the vticks are
      the inverse of the desired rates for the flows in bits per second.
  """
  def __init__(self, env, rate, vticks, debug=False):
    self.env = env
    self.rate = rate
    self.vticks = vticks
    self.auxVCs = [0.0 for i in range(len(vticks))]  # Initialize all the auxVC variables
    self.out = None
    self.n_recved = 0
    self.n_dropped = 0
    self.debug = debug
    self.store = StampedStore(env)
    self.action = env.process(self.run())  # starts the run() method as a SimPy process

  def run(self):
    while True:
      msg = (yield self.store.get())
      # Send message
      yield self.env.timeout(msg.size*8.0/self.rate)
      self.out.put(msg)

  def put(self, pkt):
    self.n_recved += 1
    now = self.env.now
    flow_id = pkt.flow_id
    # Update of auxVC for the flow. We assume that vticks is the desired bit time
    # i.e., the inverse of the desired bits per second data rate.
    # Hence we then multiply this value by the size of the packet in bits.
    self.auxVCs[flow_id] = max(now, self.auxVCs[flow_id]) + self.vticks[flow_id]*pkt.size*8.0
    # Lots of work to do here to implement the queueing discipline
    return self.store.put((self.auxVCs[flow_id], pkt))


class WFQServer(object):
  """ Models a WFQ/PGPS server. For theory and implementation see:
    Parameters
    ----------
    env : simpy.Environment
      the simulation environment
    rate : float
      the bit rate of the port
    phid : A list
      list of the phis parameters (for each possible packet flow_id). We assume a simple assignment of
      flow id to phis, i.e., flow_id = 0 corresponds to phis[0], etc...
  """
  def __init__(self, env, rate, phis, debug=False):
    self.env = env
    self.rate = rate
    self.phis = phis
    self.F_times = [0.0 for i in range(len(phis))]  # Initialize all the finish time variables
    # We keep track of the number of packets from each flow in the queue
    self.flow_queue_count = [0 for i in range(len(phis))]
    self.active_set = set()
    self.vtime = 0.0
    self.out = None
    self.n_recved = 0
    self.n_dropped = 0
    self.debug = debug
    self.store = StampedStore(env)
    self.action = env.process(self.run())  # starts the run() method as a SimPy process
    self.last_update = 0.0

  def run(self):
    while True:
      msg = (yield self.store.get())
      self.last_update = self.env.now
      flow_id = msg.flow_id
      # update information about flow items in queue
      self.flow_queue_count[flow_id] -= 1
      if self.flow_queue_count[flow_id] == 0:
        self.active_set.remove(flow_id)
      # If end of busy period, reset virtual time and reinitialize finish times.
      if len(self.active_set) == 0:
        self.vtime = 0.0
        for i in range(len(self.F_times)):
          self.F_times[i] = 0.0
      # Send message
      yield self.env.timeout(msg.size*8.0/self.rate)
      self.out.put(msg)

  def put(self, pkt):
    self.n_recved += 1
    now = self.env.now
    flow_id = pkt.flow_id
    self.flow_queue_count[flow_id] += 1
    self.active_set.add(flow_id)
    phi_sum = 0.0
    for i in self.active_set:
      phi_sum += self.phis[i]
    self.vtime += (now-self.last_update)/phi_sum
    self.F_times[flow_id] = max(self.F_times[flow_id], self.vtime) + pkt.size*8.0/self.phis[flow_id]
    # print "Flow id = {}, packet_id = {}, F_time = {}".format(flow_id, pkt.id, self.F_times[flow_id])
    self.last_update = now
    return self.store.put((self.F_times[flow_id], pkt))