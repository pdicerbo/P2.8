all: parallel

serial:
	$(MAKE) -C src
clean:
	$(MAKE) -C src clean
parallel:
	$(MAKE) -C src parallel

data:
	rm -f trajectory.xyz output.xyz energies.dat	
