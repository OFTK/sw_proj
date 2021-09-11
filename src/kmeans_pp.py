import sys, argparse
import pandas as pd
import numpy as np
import mykmeanssp

_datatype = np.float64
_output_percision = 4


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


def main():

    np.random.seed(0)

    # argument parsing
    parser = argparse.ArgumentParser(description='process kmeans++ inputs.')
    parser.add_argument('k', type=int, nargs=1, help='number of clusters')
    parser.add_argument('max_iter', type=int, nargs='?', default=300,
        help='max number of iterations')
    parser.add_argument('file_name_1', metavar='file_name_1', type=str, nargs=1, 
        help='first input file')
    parser.add_argument('file_name_2', metavar='file_name_2', type=str, nargs=1, 
        help='second input file')
    args = parser.parse_args()

    k = args.k[0]
    file_name_1 = args.file_name_1[0]
    file_name_2 = args.file_name_2[0]
    max_iter = args.max_iter

    # read files, merge DataFrames, adjust table correctly and sort by index 
    df1 = pd.read_csv(file_name_1, header=None, dtype=_datatype)
    df2 = pd.read_csv(file_name_2, header=None, dtype=_datatype)
    df = df1.merge(df2, how='inner', right_on=0, left_on=0).set_index(0, drop=True).sort_index()
    df.index = df.index.map(int)

    dim = df.shape[1]
    num_of_datapoints = df.shape[0]

    assert (k < len(df)), "ERROR: k is equal or greater than n (number of observations)"

    # run the intial centroids chioce algorithm
    initial_centroids_index = kmeans_pp(df, k)

    # print centroids correctly
    for i in range(len(initial_centroids_index)-1):
        print(initial_centroids_index[i], end=",")
    print(initial_centroids_index[len(initial_centroids_index)-1])

    # run kmeans implementation in c
    final_centroids_arr = mykmeanssp.fit(
        dim, k, max_iter,
        df.loc[initial_centroids_index].to_numpy().reshape(1, (k*dim)).tolist()[0],
        df.to_numpy().reshape(1, (num_of_datapoints*dim)).tolist()[0]
    )

    print_finalized_array(final_centroids_arr, dim)

if __name__ == "__main__":
    main()