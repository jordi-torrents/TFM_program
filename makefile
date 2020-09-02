FLAGS=-Wall -fbounds-check -g
# OPT=-O3# -freal-8-real-4
COMP=gfortran
INP=input.dat
OBJ=obj-mod/
SRC=code/
EXE=levy_prog.exe

$(EXE) : $(OBJ)def_variables.o $(OBJ)read_input.o $(OBJ)init.o $(OBJ)integrators.o $(OBJ)various.o $(OBJ)main.o $(OBJ)r1279.o $(OBJ)ran2.o $(OBJ)measurements.o $(OBJ)interaction.o
	$(COMP) $(FLAGS) $(OPT) $(OBJ)def_variables.o $(OBJ)read_input.o $(OBJ)init.o $(OBJ)integrators.o $(OBJ)various.o $(OBJ)main.o $(OBJ)r1279.o $(OBJ)ran2.o $(OBJ)measurements.o $(OBJ)interaction.o -o $(EXE)

$(OBJ)r1279.o : $(SRC)r1279/r1279.f90
	$(COMP) $(OPT) -c $< -o $@

$(OBJ)ran2.o : $(SRC)r1279/ran2.f
	$(COMP) $(OPT) -c $< -o $@

$(OBJ)def_variables.o : $(SRC)def_variables.f90
	mkdir -p $(OBJ)
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

$(OBJ)init.o : $(SRC)init.f90 $(OBJ)def_variables.o $(OBJ)integrators.o $(OBJ)read_input.o $(OBJ)various.o
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

$(OBJ)various.o : $(SRC)various.f90 $(OBJ)def_variables.o
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

$(OBJ)interaction.o : $(SRC)interaction.f90 $(OBJ)def_variables.o
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

$(OBJ)integrators.o : $(SRC)integrators.f90 $(OBJ)def_variables.o $(OBJ)interaction.o
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

$(OBJ)measurements.o : $(SRC)measurements.f90 $(OBJ)def_variables.o
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

$(OBJ)read_input.o : $(SRC)read_input.f90 $(OBJ)def_variables.o
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

$(OBJ)main.o : $(SRC)main.f90 $(OBJ)def_variables.o $(OBJ)read_input.o $(OBJ)init.o $(OBJ)integrators.o $(OBJ)various.o $(OBJ)measurements.o $(OBJ)interaction.o
	$(COMP) $(FLAGS) $(OPT) -J $(OBJ) -c $< -o $@

.PHONY: run
run:
	./$(EXE) $(INP)

.PHONY: nohup
nohup:
	nohup ./$(EXE) $(INP) &

.PHONY: clean
clean:
	rm -r $(OBJ)

.PHONY: cleanAll
cleanAll:
	rm -r $(OBJ)
	rm $(EXE)
