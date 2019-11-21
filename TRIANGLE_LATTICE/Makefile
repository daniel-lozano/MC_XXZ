## ssh dlozano@cedar.computecanada.ca
## scp <file> dlozano@cedar.computecanada.ca:<direction>
## Pass size of the 2D array with ARGS
.PHONY:plotea

plotea:Energy_magnetization.txt plot.py
	python plot.py ${ARGS}

Energy_magnetization.txt: MC_H.out
	./MC_H.out ${ARGS} ${Jperp}

MC_H.out:MC_XXZ_triangle.cpp
	g++ MC_XXZ_triangle.cpp -o MC_H.out
##	g++ -std=c++11 MC_XXZ_triangle.cpp -o MC_H.out

clean:
	rm Energy*.txt *.out
