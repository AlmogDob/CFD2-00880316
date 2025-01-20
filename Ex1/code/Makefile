CFLAGS = -Wall -Wextra -lm -fsanitize=address

solver: build_solver run_solver clean_solver  
	@echo solver done

build_solver: solver.c
	@gcc ./solver.c $(CFLAGS) -o solver

run_solver:
	@echo
	@./solver $(IN_FILE) $(OUT_DIR) 

clean_solver:
	@echo
	rm solver

debug_build_solver: solver.c
	@gcc ./solver.c $(CFLAGS) -g -o solver

##########################################################
automat: build_automat run_automat clean_automat  
	@echo automat done

build_automat: automat.c
	@gcc ./automat.c $(CFLAGS) -o automat

run_automat:
	@echo
	@./automat $(IN_FILE) $(OUT_DIR) 

clean_automat:
	@echo
	rm automat
