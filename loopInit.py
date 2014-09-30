import numpy as np

def loopInfo( numberOfReactions, input, output):
    """
    Convert the initial data in efficient data for the Cython part of the simulation

    input
        input: 2x2 matrix (Reactions | Reactants --). with the reactants and number of molecules needed for a  reaction.
        output: 2x2 matrix (Reactions | Reactants --). with the reactants and number of molecules resulting from a reaction.
        numberOfReactions: Number of reactions

    output
        updateNmatrix: input - output: 2x increase of efficiency in the molecules (N) refresh (When a reaction occurs newN = oldN - output[reaction])
        vecMol2: Vector of index of reactions to refresh when a given reaction occurs (eg: Reaction 1: 1 --> 2 + 3 //Reaction 2: 3 --> 4 // Reaction 3: 5 --6 /// Reaction 1 affects the number of molecules of reactions 1 and 2)
        vecN: Index of reactants involved in a reaction (eg: for reaction 2: 3)

    Javier Garcia-Bernardo. Oct-2012
    """


    # Update N matrix
    updateNmatrix = input - output

    # Dependende graph (Modify one reaction affects the a of the ones in the  matrix)
    affectedByReaction = input + output

    # Vector of reaction we need to refresh for each reaction (If reaction 1 affects reactions 1 2 3 and 6 --> 1 1 1 0 0 1 ...
    vecMol = np.dot(input, np.transpose(affectedByReaction))
    vecMol = np.transpose(vecMol)

    # Get the index of the reactions (0,1,2 and 5 for the example above)
    vecMol2 = list()
    for vecCount in range(numberOfReactions):
        vecMol2.append(np.nonzero(vecMol[vecCount,:]))

    # To update h quicker
    vecN = list()
    for vecNCount in range(numberOfReactions):
        vecN.append(np.nonzero(input[vecNCount,:]))

    return [updateNmatrix, vecMol2, vecN]