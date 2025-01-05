solver: build_solver run_solver clean_solver  
	@echo solver done

build_solver: solver.c
	@gcc -Wall -Wextra ./solver.c -lm -o solver

run_solver:
	@echo
	@./solver $(IN_FILE) $(OUT_DIR) 

clean_solver:
	@echo
	rm solver

debug_build_solver: solver.c
	@gcc -Wall -g -Wextra ./solver.c -lm -o solver
