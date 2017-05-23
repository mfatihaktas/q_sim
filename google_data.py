import matplotlib, numpy, pprint
# matplotlib.rcParams['pdf.fonttype'] = 42
# matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import gzip, csv
from collections import namedtuple
from rvs import *

from patch import *

"""
task events table contains the following fields:
1. timestamp
2. missing info
3. job ID
4. task index - within the job
5. machine ID
6. event type
7. user name
8. scheduling class
9. priority
10. resource request for CPU cores
11. resource request for RAM
12. resource request for local disk space
13. different-machine constraint
"""

f_to_i = {
  'timestamp': 0,
  'job_i': 2,
  'task_i': 3,
  'event': 5
}
e_to_i = {
  'schedule': 1,
  'finish': 4
}

def counter_to_furl(counter):
  part = str(counter)
  part = (5 - len(part) )*'0' + part
  # return "/home/ubuntu/task_events/part-" + part + "-of-00500.csv.gz"
  return "/cac/u01/mfa51/Desktop/google_cluster_data/task_events/part-" + part + "-of-00500.csv.gz"

def deneme():
  job_task_i__sch_finish_time_m = {}
  counter = 0
  while 1:
    furl = counter_to_furl(counter)
    try:
      with gzip.open(furl, mode="rt") as f:
        reader = csv.reader(f)
        for line in reader:
          i = line[f_to_i['job_i'] ] + '_' + line[f_to_i['task_i'] ]
          e = int(line[f_to_i['event'] ] )
          if e == e_to_i['schedule'] or e == e_to_i['finish']:
            t = float(line[f_to_i['timestamp'] ] )/10**6
            if i not in job_task_i__sch_finish_time_m:
              job_task_i__sch_finish_time_m[i] = [t]
            else:
              job_task_i__sch_finish_time_m[i].append(t)
    except (OSError, IOError) as e:
      log(WARNING, "done with the files.")
      break
    counter += 1
    if counter > 10:
      break
  with open("task_lifetime.dat", 'wt') as f:
    writer = csv.writer(f, delimiter=',')
    for job_task_i,sch_finish_time in job_task_i__sch_finish_time_m.items():
      if len(sch_finish_time) >= 2:
        sch_finish_time = [t for t in sch_finish_time if t]
        if len(sch_finish_time) == 1:
          sch_finish_time.append(0)
        # elif len(sch_finish_time) > 2:
        #   log(WARNING, "More than 2 scheduling or finish events for single task; sch_finish_time= {}".format(sch_finish_time) )
        lifetime = abs(sch_finish_time[1] - sch_finish_time[0] )
        writer.writerow([job_task_i, lifetime] )

def write_num_tasks_per_job():
  wf = open("num_tasks.dat", 'wt')
  writer = csv.writer(wf, delimiter=',')
  counter = 0
  while 1:
    print("counter= {}".format(counter) )
    ji__ti_l_m = {}
    furl = counter_to_furl(counter)
    try:
      with gzip.open(furl, mode="rt") as f:
        reader = csv.reader(f)
        for line in reader:
          ji = int(line[f_to_i['job_i'] ] )
          ti = int(line[f_to_i['task_i'] ] )
          e = int(line[f_to_i['event'] ] )
          if e == e_to_i['schedule']:
            if ji not in ji__ti_l_m:
              ji__ti_l_m[ji] = set()
            ji__ti_l_m[ji].add(ti)
      print("counter= {}, writing now...".format(counter) )
      for ji,ti_l in ji__ti_l_m.items():
        writer.writerow([ji, len(ti_l) ] )
    except (OSError, IOError) as e:
      log(WARNING, "done with the files.")
      break
    counter += 1
    if counter > 510:
      break
  wf.close()

def do_possible_merges_in_num_tasks():
  ji__num_task_m = {}
  with open("num_tasks.dat", mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      ji = int(line[0] )
      num_task = int(line[1] )
      if ji not in ji__num_task_m:
        ji__num_task_m[ji] = 0
      ji__num_task_m[ji] += num_task
  with open("num_tasks_merged.dat", mode="wt") as f:
    writer = csv.writer(f, delimiter=',')
    for ji,num_tasks in ji__num_task_m.items():
      writer.writerow([ji, num_tasks] )
  log(WARNING, "done.")

def write_jobs_w_num_task(num_task):
  ji_l = []
  with open("num_tasks_merged.dat", mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      num_task_ = int(line[1] )
      if num_task_ == num_task:
        ji_l.append(int(line[0] ) )
  print("writing, len(ji_l)= {}".format(len(ji_l) ) )
  with open("jobs_w_num_task_{}.dat".format(num_task), mode="wt") as f:
    writer = csv.writer(f, delimiter=',')
    for ji in ji_l:
      writer.writerow([ji] )
  log(WARNING, "done.")

def write_task_lifetimes(num_task):
  ji_l = []
  with open("jobs_w_num_task_{}.dat".format(num_task), mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      ji_l.append(int(line[0] ) )
  # 
  Entry = namedtuple('Entry', 'ji ti')
  entry__sch_fin_l_m = {}
  counter = 0
  while 1:
    print("counter= {}".format(counter) )
    furl = counter_to_furl(counter)
    try:
      with gzip.open(furl, mode="rt") as f:
        reader = csv.reader(f)
        for line in reader:
          ji = int(line[f_to_i['job_i'] ] )
          if ji in ji_l:
            e = int(line[f_to_i['event'] ] )
            if e == e_to_i['schedule'] or e == e_to_i['finish']:
              ti = int(line[f_to_i['task_i'] ] )
              entry = Entry(ji=ji, ti=ti)
              t = float(line[f_to_i['timestamp'] ] )/10**6
              if entry not in entry__sch_fin_l_m:
                entry__sch_fin_l_m[entry] = [0,0]
              if e == e_to_i['schedule']:
                entry__sch_fin_l_m[entry][0] = t
              elif e == e_to_i['finish']:
                entry__sch_fin_l_m[entry][1] = t
    except (OSError, IOError) as e:
      log(WARNING, "done with the files.")
      break
    counter += 1
    if counter > 510:
      break
  print("writing now...")
  with open("task_lifetimes_for_jobs_w_num_task_{}.dat".format(num_task), mode="wt") as f:
    writer = csv.writer(f, delimiter=',')
    for entry,sch_fin_tuple in entry__sch_fin_l_m.items():
      if sch_fin_tuple[0] < sch_fin_tuple[1]:
        writer.writerow([sch_fin_tuple[1]-sch_fin_tuple[0] ] )
  log(WARNING, "done.")

def plot_num_tasks_hist():
  num_tasks_l = []
  with open("num_tasks_merged.dat", mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      num_task = int(line[1] )
      # if num_task > 1000:
      #   print("num_task= {}".format(num_task) )
      # if num_task > 1 and num_task < 2000:
      num_tasks_l.append(num_task)
  
  num_task__num_job_m = {}
  for n in num_tasks_l:
    if n not in num_task__num_job_m:
      num_task__num_job_m[n] = 0
    num_task__num_job_m[n] += 1
  print("num_task__num_job_m= {}".format(pprint.pformat(num_task__num_job_m) ) )
  # plot.hist(num_tasks_l, bins=20, histtype='step')
  # plot.xlabel("Number of tasks")
  # plot.ylabel("Frequency")
  # plot.savefig("plot_num_tasks_hist.png", bbox_inches='tight')
  # plot.gcf().clear()
  log(WARNING, "done.")

def plot_task_lifetime_hist(k):
  lifetime_l = []
  with open("task_lifetimes_for_jobs_w_num_task_{}.dat".format(k), mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      lifetime_l.append(float(line[0] ) )
  
  # rv = Pareto(a=2, loc=2)
  # for i in range(1000000):
  #   lifetime_l.append(rv.gen_sample() )
  lifetime_l = numpy.sort(lifetime_l)
  print("len(lifetime_l)= {}".format(len(lifetime_l) ) )
  
  fig = plot.figure(1)
  # def_size = fig.get_size_inches()
  # fig.set_size_inches(def_size[0]*1.5, def_size[1] )
  plot.subplot(211)
  # plot.step(x_l, y_l, 'bo', label='log-linear', lw=2)
  plot.hist(lifetime_l, bins=100, histtype='step', normed=True, lw=2)
  plot.xlabel("X (s)")
  plot.ylabel("Frequency")
  plot.title(r'$k= {}$'.format(k) )
  
  x_l = lifetime_l[::-1]
  y_l = numpy.arange(lifetime_l.size)/lifetime_l.size
  plot.subplot(223)
  plot.yscale('log')
  plot.step(x_l, y_l, 'bo', label='log(tail)-X', lw=2)
  plot.xlabel("X (s)")
  plot.ylabel("Tail")
  plot.legend()
  plot.subplot(224)
  plot.xscale('log')
  plot.yscale('log')
  plot.step(x_l, y_l, 'bo', label='log(tail)-log(X)', lw=2)
  plot.xlabel("X (s)")
  plot.legend()
  
  # plot.xlabel("X")
  # plot.xlabel("Task lifetime X (s)")
  # plot.ylabel(r'$Pr\{X > x\}$')
  plot.savefig("plot_task_lifetime_hist_k_{}.png".format(k) )
  plot.gcf().clear()
  log(WARNING, "done.")

def pplot_task_lifetime_hist(k):
  lifetime_l = []
  with open("task_lifetimes_for_jobs_w_num_task_{}.dat".format(k), mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      lifetime_l.append(float(line[0] )/60/60)
  
  lifetime_l = numpy.sort(lifetime_l)
  print("len(lifetime_l)= {}".format(len(lifetime_l) ) )
  # 
  # results, edges = numpy.histogram(lifetime_l, bins=50, normed=True)
  # binWidth = edges[1] - edges[0]
  # plot.bar(edges[:-1], results*binWidth, binWidth)
  
  # plot.xlabel("Task execution time (hour)")
  # plot.ylabel("Fraction of tasks")
  # plot.title(r'Execution trace of the same job with ${}$ tasks'.format(k) )
  
  x_l = lifetime_l[::-1]
  y_l = numpy.arange(lifetime_l.size)/lifetime_l.size
  # plot.yscale('log')
  plot.step(x_l, y_l, 'bo', lw=1, linestyle=':')
  # plot.xlabel(r'$x$ (hour)')
  # plot.ylabel(r'$Pr\{X > x\}$')
  plot.xlabel(r'Task lifetime (hour)', fontsize=13)
  plot.ylabel(r'Tail distribution', fontsize=13)
  # plot.title(r'Empirical tail distribution of task exection time for jobs with {} tasks'.format(k) )
  plot.title(r'Jobs with {} tasks'.format(k), fontsize=13)
  
  fig = plot.gcf()
  def_size = fig.get_size_inches()
  print("def_size= {}".format(def_size) )
  # fig.set_size_inches(def_size[0]/1.4, def_size[1]/1.4)
  fig.set_size_inches(4, 3)
  fig.tight_layout()
  # plot.savefig("pplot_task_lifetime_hist_k_{}.pdf".format(k) )
  plot.savefig("pplot_task_lifetime_hist_k_{}.png".format(k) )
  fig.clear()

if __name__ == "__main__":
  # Uncomment with caution!
  # write_num_tasks_per_job()
  # do_possible_merges_in_num_tasks()
  
  # write_jobs_w_num_task(num_task=15)
  # write_jobs_w_num_task(num_task=400)
  # write_jobs_w_num_task(num_task=1000)
  # write_jobs_w_num_task(num_task=1050)
  
  # write_task_lifetimes(num_task=15)
  # write_task_lifetimes(num_task=400)
  # write_task_lifetimes(num_task=1000)
  # write_task_lifetimes(num_task=1050)
  
  # plot_num_tasks_hist()
  # plot_task_lifetime_hist(k=15)
  # plot_task_lifetime_hist(k=400)
  # plot_task_lifetime_hist(k=1000)
  # plot_task_lifetime_hist(k=1050)
  
  pplot_task_lifetime_hist(k=400)
  pplot_task_lifetime_hist(k=1000)
  pplot_task_lifetime_hist(k=1050)