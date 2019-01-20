
all: gcc run

gcc: main.c mymetods.c task.c
	gcc -o exefile main.c mymetods.c task.c -fopenmp -std=c11 -lm
	chmod +x exefile

run:
	./exefile

upload:
	sshpass -e scp * miet004@imm10.keldysh.ru:~/Desktop/finalprog

uploadmain:
	sshpass -e scp main.c miet004@imm10.keldysh.ru:~/Desktop/finalprog

clean:
	-rm exefile;
