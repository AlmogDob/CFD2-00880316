solver: build_solver run_solver clean_solver  
	@echo solver done

build_solver: solver.c
	gcc -g -Wall -Wextra ./solver.c -lm -o solver

run_solver:
	./solver $(IN_FILL) $(OUT_DIR) 

clean_solver:
	rm solver
