import gzip, csv
from collections import namedtuple

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
  return "/home/ubuntu/task_events/part-" + part + "-of-00500.csv.gz"

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
  with open("task_lifetime", 'wt') as f:
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
  # Entry = namedtuple('Entry', 'ji ti')
  wf = open("num_tasks", 'wt')
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
  with open("num_tasks", mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      ji = int(line[0] )
      num_task = int(line[1] )
      if ji not in ji__num_task_m:
        ji__num_task_m[ji] = 0
      ji__num_task_m[ji] += num_task
  with open("num_tasks_merged", mode="wt") as f:
    writer = csv.writer(f, delimiter=',')
    for ji,num_tasks in ji__num_task_m.items():
      writer.writerow([ji, num_tasks] )
  log(WARNING, "done.")

def plot_num_tasks_histogram():
  

if __name__ == "__main__":
  # deneme()
  # write_num_tasks_per_job()
  # do_possible_merges_in_num_tasks()
  
  
  
  
  
  
  
  
  