Get S/N in emission, abs, break features.

List from sgnaps will do.

(how to handle features that can be either? - leave status until EW is measured)

Get classifications from observers.
Machine learn the rank ordering of observability.
Features: breaks, S/N em1, S/N em2 etc in order of S/N, plus abs. Unmeasured have S/N=0.
Ignores complexity of human intuition, but we cannot really solve that.


To get S/N. We have from Chihway the S/N spectrum. Could use that directly?
No, probably not, because it doesn't really account for the strengh of the line.
Do we want an EW threshold too? Not really.... but something along those lines needs to be there. Well, it's the S/N in the line, not the cont. should be fine!
OK, so need to estimate continuum anyway - heading towards a normal fit.

Use python mpfit to fit a gaussian plus continuum (local linear slope) over some interval - get formal error out -> S/N in line.
Fix centre, fit ampl and width. Then have to margianlise down to flux error. Or can mpfit give us the error on the flux dircetly?

Saturated absorption features will not follow Gaussian line profiles. The scaling to S/N of integrated flux will be wrong - obviously.
