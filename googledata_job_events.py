import json

from google_data import *

def dump_ljn_info_m():
  ljn_info_m = {}
  for i in range(500):
  # for i in range(1, 5):
    print("i= {}".format(i) )
    furl = counter_to_furl(i, obj="job")
    try:
      with gzip.open(furl, mode="rt") as f:
        reader = csv.reader(f)
        for line in reader:
          ljn = line[jobevents_f_to_i['logical job name'] ]
          jn = line[jobevents_f_to_i['job name'] ]
          ji = line[jobevents_f_to_i['job id'] ]
          e = int(line[jobevents_f_to_i['event'] ] )
          if e == e_to_i['schedule']:
            if ljn not in ljn_info_m:
              ljn_info_m[ljn] = {'njob': 0, 'jinfo_l': [] }
            ljn_info_m[ljn]['njob'] += 1
            ljn_info_m[ljn]['jinfo_l'].append({'jn': jn, 'ji': ji} )
    except (OSError, IOError) as e:
      log(INFO, "done with the files.")
      break
  
  # json.dumps(ljn_info_m, indent=2, sort_keys=True)
  log(INFO, "dumping ljn_info_m...")
  with open('ljn_info_m.dat', 'w') as f:
    json.dump(ljn_info_m, f, indent=2, sort_keys=True)
  log(INFO, "done.")

def dump_ji_info_m():
  ji_info_m = {}
  for i in range(500):
  # for i in range(5):
    print("i= {}".format(i) )
    furl = counter_to_furl(i, obj="task")
    try:
      with gzip.open(furl, mode="rt") as f:
        reader = csv.reader(f)
        for line in reader:
          ji = line[taskevents_f_to_i['job id'] ]
          ti = line[taskevents_f_to_i['task index'] ]
          ts = float(line[taskevents_f_to_i['timestamp'] ] )/10**6
          e = int(line[taskevents_f_to_i['event'] ] )
          
          if ji not in ji_info_m:
            ji_info_m[ji] = {}
          if ti not in ji_info_m[ji]:
            ji_info_m[ji][ti] = []
          
          event_type = None
          if e == e_to_i['schedule']:
            event_type = 's'
          elif e == e_to_i['finish']:
            event_type = 'f'
          if event_type is not None:
            ji_info_m[ji][ti].append((event_type, ts) )
    except (OSError, IOError) as e:
      log(INFO, "done with the files.")
      break
  log(INFO, "dumping ji_info_m...")
  with open('ji_info_m.dat', 'w') as f:
    json.dump(ji_info_m, f, indent=2, sort_keys=True)
  log(INFO, "done.")

def dump_good_ljn_info_m():
  with open('ljn_info_m.dat') as f: 
    ljn_info_m = json.load(f)
  with open('ji_info_m.dat') as f: 
    ji_info_m = json.load(f)
  log(INFO, "loaded ljn_info_m and ji_info_m.")
  
  good_ljn_info_m = {}
  for ljn, info_m in ljn_info_m.items():
    if info_m['njob'] < 20:
      continue
    
    print("ljn= {}".format(ljn) )
    m = {}
    for jinfo_m in info_m['jinfo_l']:
      ji = jinfo_m['ji']
      try:
        m[ji] = ji_info_m[ji]
      except KeyError:
        pass
    good_ljn_info_m[ljn] = m
  log(INFO, "dumping good_ljn_info_m...")
  with open('good_ljn_info_m.dat', 'w') as f:
    json.dump(good_ljn_info_m, f, indent=2, sort_keys=True)
  log(INFO, "done.")

def dump_ljn_taskinfo_m():
  with open('good_ljn_info_m.dat') as f: 
    good_ljn_info_m = json.load(f)
  log(INFO, "loaded good_ljn_info_m.")
  
  ljn_taskinfo_m = {}
  for ljn, ji__ti__e_t_l_m_m in good_ljn_info_m.items():
    task_exectime_m = {}
    for ji, ti__e_t_l_m in ji__ti__e_t_l_m_m.items():
      for ti, e_t_l in ti__e_t_l_m.items():
        st_l, ft_l = [], []
        for e, t in e_t_l:
          if e == 's':
            st_l.append(t)
          elif e == 'f':
            ft_l.append(t)
        try:
          task_exectime_m["{}_{}".format(ji, ti) ] = max(ft_l) - min(st_l)
        except ValueError:
          pass
    ljn_taskinfo_m[ljn] = task_exectime_m
  log(INFO, "dumping ljn_taskinfo_m...")
  with open('ljn_taskinfo_m.dat', 'w') as f:
    json.dump(ljn_taskinfo_m, f, indent=2, sort_keys=True)
  log(INFO, "done.")

def plot_tasklifetime_dists():
  with open('ljn_taskinfo_m.dat') as f: 
    ljn_taskinfo_m = json.load(f)
  log(INFO, "loaded ljn_taskinfo_m.")
  
  m = {}
  for ljn, taskinfo_m in ljn_taskinfo_m.items():
    # if len(taskinfo_m) > 30*1000:
    #   m[ljn] = taskinfo_m
    lt_l = [lt for _, lt in taskinfo_m.items() ]
    try:
      if max(lt_l)/min(lt_l) > 1000: # 100
        m[ljn] = taskinfo_m
    except ValueError:
      pass
  ljn_taskinfo_m = m
  
  nrows = len(ljn_taskinfo_m) # 40
  fig, axs = plot.subplots(nrows, 1)
  figsize = (5, nrows*3)
  i = 0
  for ljn, taskinfo_m in ljn_taskinfo_m.items():
    if i == nrows:
      break
    ax = axs[i]
    plot.sca(ax)
    i += 1
    
    print(">> ljn= {}".format(ljn) )
    lt_l = numpy.sort([lt for _, lt in taskinfo_m.items() ] )
    print("len(lt_l)= {}".format(len(lt_l) ) )
    
    x_l = lt_l[::-1]
    y_l = numpy.arange(lt_l.size)/lt_l.size
    plot.step(x_l, y_l, 'bo', mew=0, linestyle=':')
    
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel(r'Task lifetime (s)', fontsize=13)
    plot.ylabel(r'Tail distribution', fontsize=13)
    # plot.ylabel(r'Fraction of tasks completed in x')
    # plot.title(r'Jobs with {} tasks'.format(k), fontsize=13)
    plot.title('ljn= {}'.format(ljn), fontsize=13)
  
  plot.subplots_adjust(wspace=1, hspace=0.5)
  fig.set_size_inches(figsize[0], figsize[1] )
  plot.savefig("plot_tasklifetime_dists.pdf", bbox_inches='tight')
  plot.gcf().clear()
  log(WARNING, "done.")

def plot_tasklifetime_dist(ljn_k_l):
  with open('ljn_taskinfo_m.dat') as f: 
    ljn_taskinfo_m = json.load(f)
  log(INFO, "loaded ljn_taskinfo_m.")
  
  for ljn, k in ljn_k_l:
    print("ljn= {}, k= {}".format(ljn, k) )
    taskinfo_m = ljn_taskinfo_m[ljn]
    lt_l = [lt for _, lt in taskinfo_m.items() ]
    add_tail_dist(plot.gca(), lt_l)
    
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel(r'Task lifetime (s)', fontsize=13)
    plot.ylabel(r'Tail distribution', fontsize=13)
    plot.title('k= {}'.format(k), fontsize=13)
    plot.savefig("plot_tasklifetime_dist_{}_k{}.png".format(ljn, k), bbox_inches='tight')
    plot.gcf().clear()
  log(WARNING, "done.")

if __name__ == "__main__":
  ## Uncomment with caution!
  # dump_ljn_info_m()
  # dump_ji_info_m()
  
  # dump_good_ljn_info_m()
  # dump_ljn_taskinfo_m()
  
  plot_tasklifetime_dists()
