import argparse
import pandas as pd
import numpy as np
from enum import IntEnum
import myspkmeanssp as spksp

_datatype = np.float64

class Goal(IntEnum):
    spk = 0
    wam = 1
    ddg = 2
    lnorm = 3
    jacobi = 4

    # In order to avoid raising a 'KeyError' while actually getting a 
    # 'ValueError'
    @staticmethod
    def from_string(s):
        try:
            return Goal[s]
        except KeyError:
            raise ValueError()
    
    # In order to print each option's name correctly in the help menu
    def __str__(self):
        return self.name

def kmeans_pp(df, k):
    z = 1

    # initialise centroids DataFrame
    centroids = pd.DataFrame(columns = df.columns).append(
        df.iloc[np.random.choice(len(df))])


    while (z < k):
        z += 1

        # calc the minimum distance for each observation
        min_dist_vector = np.array(
            [np.min(np.sum(np.square(
                centroids.to_numpy() - df.iloc[i].to_numpy()
            ), axis=1), axis=0)
            for i in range(df.shape[0])])

        # calculate probability for each observation for getting selected
        # and choose accodringly
        probability_vector = min_dist_vector / np.sum(min_dist_vector)

        centroids = centroids.append(
            df.iloc[np.random.choice(len(df), p=probability_vector)])

    return [int(x) for x in centroids.index]

def main():
    
    # Reseting the random seed
    np.random.seed(0)

    # Argument parsing
    parser = argparse.ArgumentParser(description='process spkmeans inputs.')
    parser.add_argument('k', type=int, nargs=1, help='number of clusters')
    parser.add_argument('goal', type=Goal.from_string, choices=list(Goal), 
        help='function to be used, choose from [spk, wam, ddg, lnorm, jacobi]')
    parser.add_argument('file_name', metavar='file_name', type=str, nargs=1, 
        help='The input file (.txt or .csv)')
    args = parser.parse_args()

    k = args.k[0]
    goal = args.goal
    file_name = args.file_name[0]

    # Read input file, then adjust and sort the read dataframe
    df = pd.read_csv(file_name, header=None, dtype=_datatype)

    dim = df.shape[1]
    num_of_datapoints = df.shape[0]

    assert (k < len(df)), "An Error Has Occured"

    # Performing the received subtask
    goal_return_value = spksp.perform_subtask(
        dim, k, goal,
        df.to_numpy().reshape(1, (num_of_datapoints*dim)).tolist()[0])

    if goal == Goal.spk: # It means that we need to perform 
                         # the rest of the spkmeans algorithm
        
        # If we received k = 0, we need to use the hueristic
        if k == 0:
            k = int(len(goal_return_value) / num_of_datapoints)

        # Shaping into an actual df
        df = pd.DataFrame(np.array(goal_return_value)
                                    .reshape(num_of_datapoints, k))
        
        # Run the intial centroids chioce algorithm
        initial_centroids_index = kmeans_pp(df, k)

        # print centroids correctly
        for i in range(len(initial_centroids_index)-1):
            print(initial_centroids_index[i], end=",")
        print(initial_centroids_index[len(initial_centroids_index)-1])
        
        # Performing kmeans
        goal_return_value = spksp.fit(
        k, k,
        df.loc[initial_centroids_index].to_numpy().reshape(1, (k*k)).tolist()[0],
        df.to_numpy().reshape(1, (k*num_of_datapoints)).tolist()[0]
        )

if __name__ == "__main__":
    main()