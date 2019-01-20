
all: gcc run
	
gcc: main.c mymetods.c task.c
	gcc -o exefile main.c mymetods.c task.c
	chmod +x exefile

run:
	./exefile

upload:
	sshpass -e scp main.c miet004@imm10.keldysh.ru:~/Desktop/

clean:
	-rm exefile;
