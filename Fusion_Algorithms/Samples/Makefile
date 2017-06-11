all: om_sample

clean:
	rm -rf *.exe

om_sample: om_sample.c
	gcc -o om_sample om_sample.c -std=c99 -lopenmotion -lm -D_BSD_SOURCE

exec: om_sample.exe
	./om_sample -i imu_data.csv
