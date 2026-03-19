voor het gebruiken van sisl op een windows laptop gebruik google collab. Om het script van het numerieke deel te runnen in collab moet sisl geïnstalleerd zijn en de helper_functions in de file tree staan voor collab. Het toevoegen in de file tree kan gedaan worden door in de balk links (en onder opdrachten) naar 'bestanden' te gaan. Dit opent een tree met een map ... en sample data, open deze niet. Nu kan het gedownloadde bestand van de helper_functions hierin gekopieerd worden door bovenaan in de tree te gaan naar upload symbool, gebruik dit om het helper_functions.py bestand toe te voegen. 



#Uitwerkingen van a), b) en c) 
a) 
Only the pz orbitals are relevant in Hückel calculations of C-only conjugated hydrocarbons. The other orbitals are interacting strongly with the 1s orbitals of the hydrogen orbitals and will not contribute to the chemistry, hence they are neglected. The energies alpha that are relevant  are the on-site energie of the electron of a carbon atom. The relevant interaction energie beta is the hopping energie of the electron with its nearest neighbor. 

b) 
alpha_C is the energy of the electron bounded to the carbon atom and alpha_N is the energy of the electron bounded to the nitrogen atom. The nitrogen nucleus has 1 more proton and by screening will have a larger effective charge than carbon. Thus the orbitals and valence electrons will be more tightly bounded in nitrogen than in carbon. Thus alpha_N will be more negative than alpha_N, thus alpha_C is larger than alpha_N

c)
The nitrogen atom in Pyrrole has 3 neighbours and has the same orbitals as carbon. Thus we can create a new hybrid basis (sp^2) which takes up 3 electrons from the total of 5 valence electrons. As these 3 will form sigma bonds and are neglected in the Hückel calculation, the nitrogen atom will contribute 2 electrons to the pi-system. 
