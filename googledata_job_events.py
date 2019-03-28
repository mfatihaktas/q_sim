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
      log(WARNING, "done with the files.")
      break
  
  # json.dumps(ljn_info_m, indent=2, sort_keys=True)
  log(WARNING, "dumping ljn_info_m...")
  with open('ljn_info_m.dat', 'w') as f:
    json.dump(ljn_info_m, f, indent=2, sort_keys=True)
  log(WARNING, "done.")

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
          
          if e == e_to_i['schedule']:
            event_type = 's'
          elif e == e_to_i['finish']:
            event_type = 'f'
          ji_info_m[ji][ti].append((event_type, ts) )
    except (OSError, IOError) as e:
      log(WARNING, "done with the files.")
      break
  log(WARNING, "dumping ji_info_m...")
  with open('ji_info_m.dat', 'w') as f:
    json.dump(ji_info_m, f, indent=2, sort_keys=True)
  log(WARNING, "done.")

def dump_good_ljn_info_m():
  with open('ljn_info_m.dat') as f: 
    ljn_info_m = json.load(f)
  with open('ji_info_m.dat') as f: 
    ji_info_m = json.load(f)
  
  good_ljn_info_m = {}
  for ljn, info_m in ljn_info_m.items():
    if info_m['njob'] < 20:
      continue
    
    m = {}
    for jinfo_m in info_m['jinfo_l']:
      ji = jinfo_m['ji']
      m[ji] = ji_info_m[ji]
    good_ljn_info_m[ljn] = m
  log(WARNING, "dumping good_ljn_info_m...")
  with open('good_ljn_info_m.dat', 'w') as f:
    json.dump(good_ljn_info_m, f, indent=2, sort_keys=True)
  log(WARNING, "done.")

def dump_ljn_taskinfo_m():
  with open('good_ljn_info_m.dat') as f: 
    good_ljn_info_m = json.load(f)
  
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
        task_exectime_m["{}_{}".format(ji, ti) ] = max(ft_l) - min(st_l)
    ljn_taskinfo_m[ljn] = task_exectime_m
  log(WARNING, "dumping ljn_taskinfo_m...")
  with open('ljn_taskinfo_m.dat', 'w') as f:
    json.dump(ljn_taskinfo_m, f, indent=2, sort_keys=True)
  log(WARNING, "done.")
  
if __name__ == "__main__":
  ## Uncomment with caution!
  # dump_ljn_info_m()
  dump_ji_info_m()
  
  # dump_good_ljn_info_m()
