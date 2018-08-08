# Least-Squares-Audio-Processing
<pre>
Simple audio processing techniques using least squares method
(mainly based on Ivan Selesnick's works)
(http://eeweb.poly.edu/iselesni/lecture_notes/least_squares/least_squares_SP.pdf)
</pre>
* Course project of MGE201 Operations Research I, UNIST

## Sound Files Description
* ah.wav 	: Pavarotti's High C sound (will be read only 0.10 sec for declipping_solver.m)
* beep.wav 	: Ordinal beep sound (will be read read only 0.05 sec for retrieving.m, but be played 5 times slower) 
* chu.wav	: K-POP group, "A-Pink"'s song (will be read for retrieving_solver.m)
* klaxon.wav	: Ordinal car claw sound (will be read only 0.10 sec for declipping.m)

## Matlab Codes Description
* declipping.m
	* for a sound input, make a clipped sound with user-defined cut-off value, and returns a de-clipped, and smoothed sound without using a solver (rather, uses Gram-Schmidt algorithm, QR factorization, Back Substitution)
	* -> 'declipped_and_smoothed_LS.wav', 'declipping_figure_LS' will be generated
* retrieving.m
	* for a sound input, make a missing sample with random operation, and returns a retieved sound without using a solver (rather, uses Gram-Schmidt algorithm, QR factorization, Back Substitution)
	* -> 'estimated_LS.wav', 'estimated_figure_LS' will be generated
* declipping_solver.m
	* same as 'declipping.m' EXCEPT using a built-in solver. thus, faster than model without using solver
	* -> 'declipped_and_smoothed_solver.wav', 'declipping_figure_solver' will be generated
* retrieving_solver.m
	* same as 'retrieving.m' EXCEPT using a built-in solver. thus, faster than model without using solver
	* -> 'estimated_solver.wav', 'estimated_figure_solver' will be generated
			  
## Compiling in MATLAB
* It is not recommended to just click 'RUN' button.
* Instead, you may run the code step-by-step (a 'step' is separated with '%% sign')
  you may click step and 'crtl+Enter' for compiling the step.
* 'declipping.m' and 'retrieving.m' may consume some time, (about 5 mins, 20 secs each.)
  since they use both 'get_inverse_via_GS_QR.m', and 'back_substitution.m'
  (for getting an inverse and final solution without using built-in solver)
* Check generated files (figure, sound file) to find sound input is processed well or not.

<hr>
By. Seokju Hahn / https://www.kaggle.com/ggouaeng / sjhahn11512@naver.com
