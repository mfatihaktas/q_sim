import gzip, csv

def deneme():
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
  # f = gzip.open(filename, mode="rt")
  filename = "/cac/u01/mfa51/Desktop/google_cluster_data/task_events/part-00000-of-00500.csv.gz"
  with gzip.open(filename, mode="rt") as f:
    reader = csv.reader(f)
    for line in reader:
      print("{}".format(line) )
  
if __name__ == "__main__":
  deneme()