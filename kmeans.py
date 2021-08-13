import sys

ERROR = -1

ARGV_VALID_LEN_WITHOUT_MAXITER = 2
ARGV_VALID_LEN_WITH_MAXITER = 3

def calc_vec_dist(point_a, point_b, dim):
    sum = 0
    for i in range(dim):
        sum += (point_a[i] - point_b[i])**2
    return sum

def add_dp_to_cluster(centroids, datapoint):
    min = {'value': sys.maxsize, 'centroid': 0}
    
    for centroid in centroids:
        cent_dist = calc_vec_dist(datapoint, centroid['value'], len(datapoint))
        if cent_dist < min['value']:
            min['centroid'] = centroid
            min['value'] = cent_dist
    
    min['centroid']['dp_array'].append(datapoint)

def calc_centroids_and_reset(centroids):
    dim = len(centroids[0]['value'])

    for centroid in centroids:
        sum_array = []
        for i in range(dim):
            sum_array.append(0)

        for dp in centroid['dp_array']:
            for i in range(dim):
                sum_array[i] += dp[i]
        
        for i in range(dim):
            centroid['value'][i] = sum_array[i] / len(centroid['dp_array'])

        centroid['dp_array'] = []

def kmeans(centroids, datapoints, max_iterations):
    prev_centroids = centroids
    for i in range(max_iterations):
        for dp in datapoints:
            add_dp_to_cluster(centroids, dp)

        calc_centroids_and_reset(centroids)

        n_equal = False
        for new_cent, prev_cent in zip(centroids, prev_centroids):
            if n_equal:
                break
            for n_val, p_val in zip(new_cent, prev_cent):
                if n_val != p_val:
                    n_equal = True

def read_inputs(K):
    centroids = []
    datapoints = []

    # Reading the first input separately in order to set some needed parameters such as dimension
    first_line = input().split(',')
    dim = len(first_line)

    new_point = []
    for val in first_line:
        try:
            new_point.append(float(val))
        except(ValueError):
            print(f"first datapoint's value wasn't able to be parsed as an int, exiting...")
            return ERROR

    centroids.append({'value': new_point.copy(), 'dp_array': [new_point]})
    datapoints.append(new_point)

    for line in sys.stdin:
        splitted_line = line.split(',')

        if len(splitted_line) != dim:
            print(f"All datapoints should be of dimension {dim}, received a bad input of dim {len(splitted_line)}, exiting...")
            return ERROR
        
        new_point = []
        for val in splitted_line:
            try:
                new_point.append(float(val))
            except(ValueError):
                print(f"some datapoints' value wasn't able to be parsed as an int, exiting...")
                return ERROR

        datapoints.append(new_point)

        if len(datapoints) <= K:
            centroids.append({'value': new_point.copy(), 'dp_array': [new_point]})
        else:
            add_dp_to_cluster(centroids, new_point)

    return centroids, datapoints

def parse_args():
    if sys.argv.__len__() == ARGV_VALID_LEN_WITHOUT_MAXITER:
        max_iter = 200
    elif sys.argv.__len__() == ARGV_VALID_LEN_WITH_MAXITER:
        max_iter = sys.argv.__getitem__(2)
    else:
        print(f"got {sys.argv.__len__()} arguments, expecting either 2 or 3")
        return ERROR
    
    arguments = {}
    
    try:
        arguments['K'] = int(sys.argv.__getitem__(1))
        arguments['max_iter'] = int(max_iter)

    except(ValueError):
        print(f"K ({sys.argv.__getitem__(1)}) and max_iter "
        f"({max_iter}) both should be integers")
        return ERROR
    
    if arguments['K'] < 0 or arguments['max_iter'] < 0:
        print(f"K ({arguments['K']}) and max_iter ({arguments['max_iter']})"
        " must be set to integers larger than zero")
        return ERROR
        
    return arguments

def print_centroids(centroids):
    for centroid in centroids:
        print()
        for value in centroid['value']:
            print(round(value,4), end=", ")

def main():
    args = parse_args()

    if args == ERROR:
        return
    
    centroids, datapoints = read_inputs(args['K'])
    
    kmeans(centroids, datapoints, args['max_iter'])

    print_centroids(centroids)

if __name__ == "__main__":
    main()