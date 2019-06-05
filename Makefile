all: spec_find

spec_find: spec_find.o k_input.o es_input.o find_close.o echelon.o  
	g++ -o $@ $^
	g++ -c $^ -o $@

clean:
	rm -rf generator *.o

.PHONY: all clean
