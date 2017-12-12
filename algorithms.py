import math

############################# SELECTION SORT ################################

def selection_sort(collection):
    comparisons = 0
    length = len(collection)
    for i in range(length):
        least = i
        for k in range(i + 1, length):
            comparisons += 1
            if collection[k] < collection[least]:
                least = k
        collection[least], collection[i] = (
            collection[i], collection[least]
        )
    return comparisons, collection

############################## INSERTION SORT ##############################

def insertion_sort(collection):
    comparisons = 0
    for i in range(1, len(collection)):
        tmp = collection[i]
        k = i
        while k > 0 and tmp < collection[k - 1]:
            comparisons += 1;
            collection[k] = collection[k - 1]
            k -= 1
        comparisons += 1;
        collection[k] = tmp
    return comparisons, collection

############################# BUBBLE SORT ##################################

def bubble_sort(collection):
    comparisons = 0
    while True:
        swapped = False
        for i in range(1, len(collection)):
            comparisons += 1
            if collection[i-1] > collection[i]:
                collection[i-1], collection[i] = collection[i], collection[i-1]
                swapped = True
        if not swapped:
            break
    return comparisons, collection

############################# SHELL SORT ##################################

# Using Shell Sequence
def shell_sort(arr):
    comparisons = 0
    n = len(arr)
    gap = int(n/2)
    
    while gap > 0:
 
        for i in range(gap,n):
            temp = arr[i]
            
            j = i
            while  j >= gap and arr[j-gap] > temp:
                comparisons += 1
                arr[j] = arr[j-gap]
                j -= gap
            
            arr[j] = temp
            comparisons += 1
        gap = int(gap/2)
    
    return comparisons

# Using Ciura's Gap Sequence
def shell_sort_v1(array):
    comparisons = 0
    gaps = [701, 301, 132, 57, 23, 10, 4, 1]
    
    for gap in gaps:
        
        for i in range(gap, len(array)):
            tmp = array[i]
            
            j = i
            while j >= gap and tmp < array[j - gap]:
                comparisons += 1
                array[j] = array[j - gap]
                j -= gap
                
            array[j] = tmp
            comparisons += 1
    
    return comparisons, array

######################### BUCKET SORT #################################

def bucket_sort(arr):
    # get hash codes
    comparisons = 0
    code = hashing(arr)
    buckets = [list() for _ in range(code[1])]
    # distribute data into buckets
    for i in arr:
        x = re_hashing(i, code)
        buck = buckets[x]
        buck.append(i)

    # Sort each bucket

    for bucket in buckets:
        comparisons += insertion_sort(bucket)[0]

    ndx = 0
    # merge the buckets
    for b in range(len(buckets)):
        for v in buckets[b]:
            arr[ndx] = v
            ndx += 1
    
    return comparisons, buckets
 
def hashing(arr):
    m = arr[0]
    for i in range(1, len(arr)):
        if (m < arr[i]):
            m = arr[i]
    result = [m, int(math.sqrt(len(arr)))]
    return result

def re_hashing(i, code):
    return int(i / code[0] * (code[1] - 1))

######################## RADIX SORT ###############################

# A function to do counting sort of arr[] according to
# the digit represented by exp.
def counting_sort(arr, exp1):
    accesses = 0
    n = len(arr)
 
    # The output array elements that will have sorted arr
    output = [0] * (n)
 
    # initialize count array as 0
    count = [0] * (10)
 
    # Store count of occurrences in count[]
    for i in range(0, n):
        accesses += 1
        index = int((arr[i]/exp1))
        count[ int((index)%10) ] += 1

    # Change count[i] so that count[i] now contains actual
    #  position of this digit in output array
    for i in range(1,10):
        count[i] += count[i-1]
 
    # Build the output array
    i = n-1
    while i>=0:
        accesses += 1
        index = (arr[i]/exp1)
        output[ count[ int((index)%10) ] - 1] = arr[i]
        count[ int((index)%10) ] -= 1
        i -= 1
 
    # Copying the output array to arr[],
    # so that arr now contains sorted numbers
    i = 0
    for i in range(0,len(arr)):
        arr[i] = output[i]

    return accesses

# Radix Sort Using Counting sort as Subroutine
def radix_sort_counting(arr):
    accesses = 0
    max1 = max(arr)
    exp = 1
    while int(max1/exp) > 0:
        accesses += counting_sort(arr,exp)
        exp *= 10
    return accesses

def radix_sort_insertion(arr):
    accesses = 0
    max1 = max(arr)
    exp = 1
    while int(max1/exp) > 0:
        accesses += insertion_sort_exp(arr, exp)[0]
        exp *= 10
    return accesses

def insertion_sort_exp(collection, exp):
    comparisons = 0
    for i in range(1, len(collection)):
        tmp = collection[i]
        k = i
        while k > 0 and int(tmp/exp) % 10 < int(collection[k-1]/exp) % 10:
            comparisons += 1;
            collection[k] = collection[k - 1]
            k -= 1
        comparisons += 1;
        collection[k] = tmp
    return comparisons, collection

########################### MERGE SORT #################################

# Merges two subarrays of arr[].
# First subarray is arr[l..m]
# Second subarray is arr[m+1..r]
def merge(arr, l, m, r):
    comparison = 0
    n1 = m - l + 1
    n2 = r- m
    
    # create temp arrays
    L = [0] * (n1)
    R = [0] * (n2)
 
    # Copy data to temp arrays L[] and R[]
    for i in range(0 , n1):
        L[i] = arr[l + i]
 
    for j in range(0 , n2):
        R[j] = arr[m + 1 + j]
 
    # Merge the temp arrays back into arr[l..r]
    i = 0     # Initial index of first subarray
    j = 0     # Initial index of second subarray
    k = l     # Initial index of merged subarray
 
    while i < n1 and j < n2:
        if L[i] <= R[j]:
            arr[k] = L[i]
            i += 1
        else:
            arr[k] = R[j]
            j += 1
        comparison += 1
        k += 1
 
    # Copy the remaining elements of L[], if there
    # are any
    while i < n1:
        arr[k] = L[i]
        i += 1
        k += 1
 
    # Copy the remaining elements of R[], if there
    # are any
    while j < n2:
        arr[k] = R[j]
        j += 1
        k += 1
    
    return comparison

# l is for left index and r is right index of the
# sub-array of arr to be sorted
def merge_sort(arr,l,r):
    comparison = 0;
    if l < r:
 
        # Same as (l+r)/2, but avoids overflow for
        # large l and h
        m = int((l+(r-1))/2)
 
        # Sort first and second halves
        comparison += merge_sort(arr, l, m)
        comparison += merge_sort(arr, m+1, r)
        comparison += merge(arr, l, m, r)
    return comparison

########################### QUICK SORT #############################

# This function takes last element as pivot, places
# the pivot element at its correct position in sorted
# array, and places all smaller (smaller than pivot)
# to left of pivot and all greater elements to right
# of pivot
def partition(arr,low,high):
    comparison = 0
    i = ( low-1 )         # index of smaller element
    pivot = arr[high]     # pivot
    
    for j in range(low , high):
        comparison += 1
        # If current element is smaller than or
        # equal to pivot
        if arr[j] <= pivot:
            # increment index of smaller element
            i = i+1
            arr[i],arr[j] = arr[j],arr[i]
        
    arr[i+1],arr[high] = arr[high],arr[i+1]
    return ( i+1 ), comparison

# The main function that implements QuickSort
# arr[] --> Array to be sorted,
# low  --> Starting index,
# high  --> Ending index

# Function to do Quick sort
def quick_sort(arr,low,high):
    comparison = 0

    if low < high:

        # pi is partitioning index, arr[p] is now
        # at right place
        pi, comparison = partition(arr,low,high)

        # Separately sort elements before
        # partition and after partition
        comparison += quick_sort(arr, low, pi-1)
        comparison += quick_sort(arr, pi+1, high)
        
    return comparison