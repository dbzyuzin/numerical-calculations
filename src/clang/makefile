
all: uploadmain compile run

compile: main.c mymetods.c task.c
	sshpass -e ssh miet004@imm10.keldysh.ru \
	'mpicc -o ~/Desktop/finalprog/exefile ~/Desktop/finalprog/main.c ~/Desktop/finalprog/mymetods.c ~/Desktop/finalprog/task.c -fopenmp -std=c11 -lm'

run:
	sshpass -e ssh miet004@imm10.keldysh.ru 'mpirun -np 4 ~/Desktop/finalprog/exefile'

upload:
	sshpass -e scp * miet004@imm10.keldysh.ru:~/Desktop/finalprog

uploadmain:
	sshpass -e scp main.c miet004@imm10.keldysh.ru:~/Desktop/finalprog

clean:
	-rm exefile;
