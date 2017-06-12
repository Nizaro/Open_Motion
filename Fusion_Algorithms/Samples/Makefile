all: om_sample

clean:
	rm -rf *.exe

util_.o: utils_.c utils_.h
	gcc -c  utils_.c -std=c99 -lopenmotion -lm -D_BSD_SOURCE

om_sample: om_sample.c util_.o
	gcc -o om_sample om_sample.c utils_.o -std=c99 -lopenmotion -lm -D_BSD_SOURCE

exec: om_sample.exe
	./om_sample -i imu_data.csv
