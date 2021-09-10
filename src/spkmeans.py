import argparse
import pandas as pd
import numpy as np
import myspkmeanssp as spksp

_datatype = np.float64
_output_percision = 4

class Goal(Enum):
    spk = 1
    wam = 2
    ddg = 3
    lnorm = 4
    jacobi = 5

    # In order to avoid raising a 'KeyError' while actually getting a 'ValueError'
    @staticmethod
    def from_string(s):
        try:
            return Goal[s]
        except KeyError:
            raise ValueError()

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

def print_finalized_array(arr, dim):
    ret = ""
    for i in range(len(arr)):
        if (i % dim) == (dim - 1):
            if i == (len(arr)-1): end_char = ""
            else: end_char = "\n"
        else: end_char = ","

        print(round(arr[i], _output_percision), end=end_char)

def print_according_to_goal(value, dim, goal):
    if goal == 1:
        print_finalized_array(value, dim)
    elif goal == 2:
        pass
    elif goal == 3:
        pass
    elif goal == 4:
        pass
    elif goal == 5:
        pass

def main():

    # Argument parsing
    parser = argparse.ArgumentParser(description='process spkmeans inputs.')
    parser.add_argument('k', type=int, nargs=1, help='number of clusters')
    parser.add_argument('goal', type=Goal.from_string, choices=list(Goal),
        help='function to be used, choose from [spk, wam, ddg, lnorm, jacobi]')
    parser.add_argument('file_name', metavar='file_name', type=str, nargs=1, 
        help='The input file (.txt or .csv)')
    args = parser.parse_args()

    k = args.k
    goal = args.goal
    file_name = args.file_name

    # Read input file, then adjust and sort the read dataframe
    df = pd.read_csv(file_name, header=None, dtype=_datatype).set_index(0).sort_index()
    df.index = df.index.map(int)

    dim = df.shape[1]
    num_of_datapoints = df.shape[0]

    assert (k < len(df)), "ERROR: k is equal or greater than n (number of observations)"

    if goal == 1: # It means that we need to perform the spkmeans algorithm

        # Run the intial centroids chioce algorithm
        initial_centroids_index = kmeans_pp(df, k)

        goal_return_value = spksp.spkmeans_fit(
        dim, k,
        df.loc[initial_centroids_index].to_numpy().reshape(1, (k*dim)).tolist()[0],
        df.to_numpy().reshape(1, (num_of_datapoints*dim)).tolist()[0]
        )
        
    else:

        goal_return_value = spksp.perform_subtask(
        dim, k, goal,
        df.to_numpy().reshape(1, (num_of_datapoints*dim)).tolist()[0]
        )

    print_according_to_goal(goal_return_value, dim, goal)

if __name__ == "__main__":
    main()