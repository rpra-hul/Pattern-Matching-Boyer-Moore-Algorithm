import sys


def read_file_to_string(file_path):
    """Read the contents of a text file and return it as a string.

    Helper function to read the contents of the text and pattern files and
    return them as strings to be used in the Boyer-Moore pattern matching algorithm.

    Args:
        file_path (str): The path to the text or pattern file.

    Returns:
        file_contents (str): The contents of the file as a string.
    """
    # File open in read mode
    with open(file_path, "r") as f:
        # Read the contents of the file into a string
        file_contents = f.read()
    return file_contents


def write_result_to_file(file_path, arr):
    """Write the results of the pattern matching to a text file.

    Helper function to write the starting indices of the occurences of the pattern
    to a specified output text file. The indices are 1-based indexed.

    Args:
        file_path (str): The path to the output text file.
        arr (list): The list of starting indices of the occurences of the pattern.
    """
    f = open(file_path, "w")
    for i in arr:
        f.write(str(i + 1) + "\n")
    f.close()
    return


def compute_Z_suffix_values(pat):
    """Compute the Z-suffix array of the pattern.
    
    A pre-processing step to compute the length of the longest substring 
    that matches the suffix of the pattern at each index.

    Args:
        pat (str): The pattern to search for in the text.
    
    Returns:
        Z_suffix (list):
            A list storing the length of the longest substring that matches 
            the suffix of the pattern at each index.
    """
    n = len(pat)
    # Initialise an array of Z values of size n
    Z_values = [0] * n

    # Reverse the pattern
    reverse_pat= ""
    for char in pat:
        reverse_pat = char + reverse_pat

    # Initialise left, right indices of the Z-box and an iterative variable, i
    l, r, i = 0, 0, 0
    for k in range(1, n):
        # Base Case: k = 1 or First Case: k > r, Perform explicit character comparison
        if k > r:
            i = k
            # Initialise a count variable to keep track of the number of matching characters
            count = 0
            # As long as we don't reach the end of the pattern and the characters match,
            while i < n and reverse_pat[i] == reverse_pat[i - k]:
                # Increment the count and i
                i += 1
                count += 1
            # Store the count in the Z_values array
            Z_values[k] = count
            # Update the left and right indices of the Z-box
            l = k
            r = k - 1 + count
        
        # Case 2: k <= r, k is within the Z-box
        elif k <= r:
            # Calculate the remaining characters in the Z-box
            rem = r - k + 1
            k_prime = k - l

            # Case 2a: Z[k'] < rem, Z'-box is within the remaining Z-box
            if Z_values[k_prime] < rem:
                Z_values[k] = Z_values[k_prime]

            # Case 2b: Z[k'] == rem, Perform explicit character comparison
            # starting from the end of the Z-box
            elif Z_values[k_prime] == rem:
                # Set i to the right index of the Z-box
                i = r + 1
                count = 0
                # As long as we don't reach the end of the pattern and the characters match,
                while i < n and reverse_pat[i] == reverse_pat[i - k]:
                    # Increment the count and i
                    i += 1
                    count += 1
                # Store the count in the Z_values array
                Z_values[k] = count + rem
                # Update the left and right indices of the Z-box
                l = k
                r = i - 1

            # Case 2c: Z[k'] > rem, Z'-box extends beyond the remaining Z-box 
            else:
                # Set Z[k] to the remaining characters in the Z-box
                Z_values[k] = rem

    # Reverse the Z_values array to get the Z-suffix array
    Z_suffix = []
    for i in range(len(Z_values) - 1, -1, -1):
        Z_suffix.append(Z_values[i])
    
    return Z_suffix


def compute_good_suffix(Z_suffix):
    """Compute the good suffix array based on the Z-suffix array.
    
    A pre-processing step to compute the endpoint indices of all occurrences
    of the substrings that exactly matches the suffix of the pattern.

    Args:
        Z_suffix (list[int]): The Z-suffix array of the pattern.

    Returns:
        good_suffix (list):
            A list storing the endpoint indices of all occurences of 
            the substrings that exactly matches the suffix of the pattern.
    """
    n = len(Z_suffix)
    # Initialise the good_suffix array
    good_suffix = [-1] * (n + 1)

    # Iterate throught the Z_suffix array
    for p in range(n):
        j = n - Z_suffix[p]
        # If the suffix was not present in the pattern
        if good_suffix[j] == -1:
            # Create a new array and store the endpoint index
            good_suffix[j] = [p]
        else:
            # Append the endpoint index to the existing array
            good_suffix[j].append(p)

    return good_suffix


def compute_matched_prefix(Z_suffix, n):
    """Compute the matched prefix array based on the Z-suffix array.
    
    A pre-processing step to compute the length of the largest suffix
    of the good suffix that matches the prefix of the pattern for each
    index.

    Args:
        Z_suffix (list[int]): The Z-suffix array of the pattern.
        n (int): The length of the pattern.

    Returns:
        matched_prefix (list[int]):
            A list storing the length of the largest suffix of a good suffix
            that matches the prefix of the pattern for each index.
    """
    # Initialise the matched_prefix array
    matched_prefix = [0] * (n + 1)
    # Initialise an iterative variable, j
    j = 0
    # Iterate backwards over the matched prefix array 
    for k in range(n - 1, 0, -1):
        # If the curr index is a suffix of the pattern
        if k + Z_suffix[j] == n:
            # Store the length of the matched prefix
            matched_prefix[k] = Z_suffix[j]
            j += 1
        else:
            # Copy the matched prefix value from the previous iteration
            matched_prefix[k] = matched_prefix[k + 1]
            j += 1
    # Set the length of the pattern as the matched prefix for the first index
    matched_prefix[0] = n

    return matched_prefix


def boyer_moore(txt, pat):
    """Boyer-Moore pattern matching algorithm.
    
    Performs pattern matching on the text based on a given pattern using 
    an edited Boyer-Moore algorithm with a stricter good suffix rule.

    Args:
        txt (str): The text to search for the pattern.
        pat (str): The pattern to search for in the text.
    
    Returns:
        Prints the starting indices of the occurences of the pattern in the text 
        to an output text file named "output_q1.txt". (1-based indexing)
    """
    n = len(txt)
    m = len(pat)

    # Pre-processing phase
    z_suffix = compute_Z_suffix_values(pat)
    good_suffix = compute_good_suffix(z_suffix)
    matched_prefix = compute_matched_prefix(z_suffix, m)

    # Searching phase
    result = []
    start, stop = -1, -1
    # Cumulative shift value of the pattern with respect to the text
    shift = 0
    while shift <= n - m:
        # Set i to the length of the pattern
        k = m - 1
        # Perform explicit comparison of the pattern and the text from right to left
        while k >= 0 and pat[k] == txt[k + shift]:
            # Skip over the region that is guaranteed to match in the next iteration
            if stop != -1 and (k - 1) == stop:
                k = start - 1
                start, stop = -1, -1
            else:
                # Decrement i as long as the characters match
                k -= 1

        # An occurence of the pattern is found in the text if i < 0
        if k < 0:
            # Append the starting index of the pattern in the text to the results list
            result.append(shift)
            # Shift the pattern to the right by the length of the pattern - matched_prefix[1]
            shift += m - matched_prefix[1]
        else:
            # Determine the shift value based on the good suffix rule
            # If there is no occurence of the suffix in the pattern
            if good_suffix[k + 1] == -1:
                gs_shift = m - matched_prefix[k + 1] 
                # Update start and stop values
                start = 1
                stop = matched_prefix[k + 1]
            else:
                # If there is an occurence of the suffix in the pattern
                # Store the bad character in the text
                bad_char = txt[k + shift]
                # For each occurence of the suffix in the pattern
                pos = 0
                for i in range(len(good_suffix[k + 1])):
                    # Get the starting index of the good suffix in the pattern
                    gs_value = good_suffix[k + 1][i]
                    # Get the length of the good suffix
                    z_suffix_value = z_suffix[gs_value]
                    pos = i

                    # If the position of the char to compare is to the right of the bad character
                    if (gs_value - z_suffix_value) >= k:
                        # Break out of the loop early
                        break

                    # If the character preceding the good suffix in the pattern is equal to the bad character
                    if (
                        gs_value - z_suffix_value >= 0
                        and pat[gs_value - z_suffix_value] == bad_char
                    ):
                        # Set the shift value and break out of the loop
                        gs_shift = m - gs_value - 1
                        break
                    else:
                        # Update the shift value at every iteration
                        gs_shift = m - gs_value - 1

                # Update the start and stop values
                start = good_suffix[k + 1][pos] - m + k + 1
                stop = good_suffix[k + 1][pos] 
            # Update the shift value
            shift += gs_shift

    return result


if __name__ == "__main__":
    # Retrieve the file paths from the command line arguments
    _, txt_file, pat_file = sys.argv

    # Read and store the contents of the text and pattern files
    txt_file = read_file_to_string(txt_file)
    pat_file = read_file_to_string(pat_file)

    # Perform pattern matching using the Boyer-Moore algorithm
    result = boyer_moore(txt_file, pat_file)

    # Print the starting indices of the pattern in the text to an output file
    write_result_to_file("output_q1.txt", result)
