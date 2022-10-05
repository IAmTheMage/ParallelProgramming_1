import subprocess
import os
import numpy as np
import csv

number_of_executions = int(input('Digite a quantidade de vezes que deseja fazer o teste: '))
size = int(input('Digite o tamanho do vetor: '))

os.system('mpicc -O3 -o mpi_vector_sum mpi_vector_sum.c')

executabl = 'mpiexec -n 1 ./mpi_vector_sum '
executabl_2 = 'mpiexec -n 2 ./mpi_vector_sum '
executabl_4 = 'mpiexec -n 4 ./mpi_vector_sum '
executabl_8 = 'mpiexec -n 8 ./mpi_vector_sum '
executabl_16 = 'mpiexec -n 16 ./mpi_vector_sum '
executabl += str(size)
executabl_2 += str(size)
executabl_4 += str(size)
executabl_8 += str(size)
executabl_16 += str(size)

csv_results_medians = ['Média']
csv_results_std = ['Desvio padrão']
csv_headers = ['Numero de processos','1', '2', '4', '8', '16']

zer = np.zeros(number_of_executions, dtype=np.float32)
for i in range(number_of_executions):
  print("Interação: {}".format(i))
  result = subprocess.run(executabl.split(' '), stdout=subprocess.PIPE)
  result_str = result.stdout.decode('utf-8')
  result_list = result_str.split(' ')
  result_time = float(result_list[2])
  zer[i] = result_time

print('Median of time is: {}'.format(np.median(zer)))
print('Standart deviation is: {}'.format(np.std(zer)))
csv_results_medians.append(str(np.median(zer)))
csv_results_std.append(np.std(zer))

for i in range(number_of_executions):
  print("Interação: {}".format(i))
  result = subprocess.run(executabl_2.split(' '), stdout=subprocess.PIPE)
  result_str = result.stdout.decode('utf-8')
  result_list = result_str.split(' ')
  result_time = float(result_list[2])
  zer[i] = result_time

print('Median of time is: {}'.format(np.median(zer)))
print('Standart deviation is: {}'.format(np.std(zer)))
csv_results_medians.append(str(np.median(zer)))
csv_results_std.append(np.std(zer))

for i in range(number_of_executions):
  print("Interação: {}".format(i))
  result = subprocess.run(executabl_4.split(' '), stdout=subprocess.PIPE)
  result_str = result.stdout.decode('utf-8')
  result_list = result_str.split(' ')
  result_time = float(result_list[2])
  zer[i] = result_time

print('Median of time is: {}'.format(np.median(zer)))
print('Standart deviation is: {}'.format(np.std(zer)))
csv_results_medians.append(str(np.median(zer)))
csv_results_std.append(np.std(zer))

for i in range(number_of_executions):
  print("Interação: {}".format(i))
  result = subprocess.run(executabl_8.split(' '), stdout=subprocess.PIPE)
  result_str = result.stdout.decode('utf-8')
  result_list = result_str.split(' ')
  result_time = float(result_list[2])
  zer[i] = result_time

print('Median of time is: {}'.format(np.median(zer)))
print('Standart deviation is: {}'.format(np.std(zer)))
csv_results_medians.append(str(np.median(zer)))
csv_results_std.append(np.std(zer))

for i in range(number_of_executions):
  print("Interação: {}".format(i))
  result = subprocess.run(executabl_16.split(' '), stdout=subprocess.PIPE)
  result_str = result.stdout.decode('utf-8')
  result_list = result_str.split(' ')
  result_time = float(result_list[2])
  zer[i] = result_time

print('Median of time is: {}'.format(np.median(zer)))
print('Standart deviation is: {}'.format(np.std(zer)))
csv_results_medians.append(str(np.median(zer)))
csv_results_std.append(np.std(zer))

f = open('results.csv', 'w')
writer = csv.writer(f)

writer.writerow(csv_headers)
writer.writerow(csv_results_medians)
writer.writerow(csv_results_std)
