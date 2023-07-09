# [PSI+] model from Derdowski et al., 2010, Science.
# Model utilizes Gillespie's algorithm for stochastic simulation 
# of coupled reactions in the yeast cell (see documentation for further details).
# Yury A. Barbitoff, 2015

# Loading dependencies
source('./technical.r')
source('./rSim.r')

# Provide info about simulation
cat('\n', '\n\nHello, dear User! :D', 
    '\nI`m Slowbro-code, the New, Fast, but not fancy, version of the',
    '\n[PSI+] prion model!', 'After numerous attempts to make Me look',
    '\nbetter, the Master has given up.', 'Nevertheless, he still',
    '\nmanaged to improve Me - for example, now I work up to 2.5 times',
    '\nfaster!', '\nNow let us begin our simulations! To remind you:',
    '\noutput files of the simulation are stored in your R Working Directory,',
    '\nthat is', getwd())
      


