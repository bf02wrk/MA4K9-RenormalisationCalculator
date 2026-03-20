# MA4K9-RenormalisationCalculator
A simple Python script that calculates the renormalised amplitudes of a family of simple $\phi^4$ Feynman diagrams.

This is a simple Python program I wrote to handle some of the calculations for my MA4K9 dissertation on the Connes-Kreimer theory of renormalisation. Its purpose is to calculate the renormalised amplitudes of massless $\phi^4$ Feynman graphs in the family >O<, >OO<, >OOO<, etc. 

The function calculate_renormalised_amplitude(N) does this by taking the loop number $N$ of the graph as an input, calculating the appropriate number of terms of the series expansions of $\Gamma(z)$ and $(1-2z)^{-1}$, multiplying them out, discarding unneeded terms and then applying the main theorem of Connes-Kreimer theorem to find an expression for the renormalised part of the Feynman integral in terms of $z$, $\lambda$, $L=\log\left( \frac{p^2}{4\pi \mu^2} \right)$ and the constants $\gamma$ and $\zeta(n)$, $n\geq 2$.

The function evaluate(expression, log_val, coupling_constant, loop_num) substitutes into the input expression the values $z=0$, $L=$log_val, $\lambda=$coupling_constant and uses $N=$loop_num to replace the zeta values with their numerical values. This gives us a single complex value for each loop number.

By default, the program returns the first three amplitudes and the time taken to calculate each one.

If you have any queries please don't hesitate to contact Bee Fox at foxbee02@gmail.com.
