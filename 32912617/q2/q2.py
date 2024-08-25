import sys

# Defining global variables
unknown_char = "!"
alphabet_start = 37
alphabet_length = 126 - 37 + 1


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
    m = len(pat)
    global unknown_char
    # Initialise an array of Z values of size m
    Z_values = [0] * m
    # Reverse the pattern
    reverse_pat = pat[::-1]
    # print(reverse_pat)
    # Initialise left, right indices of the Z-box and an iterative variable, i
    l, r = -1, -1
    for k in range(1, m):
        # Base Case: k = 1 or First Case: k > r, Perform explicit character comparison
        if k > r:
            unknown_char_occ = 0
            i = k
            # Initialise a count variable to keep track of the number of matching characters
            count = 0
            # As long as we dont reach the end of the pattern and the characters match,
            while i < m and (
                reverse_pat[i] == reverse_pat[i - k]
                or reverse_pat[i] == unknown_char
                or reverse_pat[i - k] == unknown_char
            ):

                if reverse_pat[i - k] == unknown_char:
                    unknown_char_occ += 1

                # Increment the count and i
                i += 1
                count += 1

            # Store the count in the Z_values array
            Z_values[k] = count
            # Update the left and right indices of the Z-box
            l = k
            r = k - 1 + count - unknown_char_occ
        # Case 2: k <= r, k is within the Z-box
        elif k <= r:
            # Calculate the remaining characters in the Z-box
            rem = r - k + 1
            k_prime = k - l
            # Case 2a: Z[k'] < rem, Z'-box is within the remaining Z-box
            if Z_values[k_prime] < rem:
                # Copy the Z[k'] value to Z[k]
                Z_values[k] = Z_values[k_prime]
            # Case 2b: Z[k'] == rem, Perform explicit character comparison
            # starting from the end of the Z-box
            elif Z_values[k_prime] == rem:
                i = r
                count = 0
                while i < m and (
                    reverse_pat[i] == reverse_pat[i - k]
                    or reverse_pat[i] == unknown_char
                    or reverse_pat[i - k] == unknown_char
                ):
                    i += 1
                    count += 1
                Z_values[k] = count
                l = k
                r = k - 1 + count
            # Case 2c: Z[k'] > rem, Z'-box extends beyond the remaining Z-box
            else:
                # Set Z[k] to the remaining characters in the Z-box
                Z_values[k] = rem
    # Reverse the Z-values to obtain the Z-suffix array
    Z_suffix = Z_values[::-1]
    return Z_suffix


def compute_good_suffix(Z_suffix):
    """Compute the good suffix array based on the Z-suffix array.
    
    A pre-processing step to compute the endpoint indices of the rightmost 
    occurences of the substrings that exactly matches the suffix of the pattern.

    Args:
        Z_suffix (list): The Z-suffix array of the pattern.

    Returns:
        good_suffix (list):
            A list storing the endpoint indices of the rightmost occurences of 
            the substrings that exactly matches the suffix of the pattern.
    """
    n = len(Z_suffix)
    # Initialise the good suffix array
    good_suffix = [-1] * (n + 1)

    #Iterate through the Z-suffix array
    for p in range(n):
        # Compute the endpoint index of the rightmost occurence of the substring
        j = n - Z_suffix[p]
        good_suffix[j] = p

    return good_suffix


def compute_matched_prefix(Z_suffix, n):
    """Compute the matched prefix array based on the Z-suffix array.
    
    A pre-processing step to compute the length of the largest suffix
    of the good suffix that matches the prefix of the pattern for each
    index.

    Args:
        Z_suffix (list): The Z-suffix array of the pattern.
        n (int): The length of the pattern.

    Returns:
        matched_prefix (list):
            A list storing the length of the largest suffix of a good suffix
            that matches the prefix of the pattern for each index.
    """
    # Initialise the matched prefix array
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


def compute_bad_char_table(pat):
    """Computes the extended bad character table of the pattern.
    
    A pre-processing step to compute the rightmost occurence of each character 
    to the left of the mismatched character in the pattern. Stored as a 2D 
    matrix of size O(m * alphabet_length) where m is the length of the pattern.

    Args:
        pat (str): The pattern to search for in the text.
    
    Returns:
        bad_char_table (list): A 2D matrix storing the index of the rightmost 
        occurence of each character to the left of the mismatched character.
    """
    # Initialise the global variables
    global unknown_char, alphabet_start, alphabet_length
    m = len(pat)
    # Create a 2D matrix to store the extended bad character table
    bad_char_table = [[-1] * (alphabet_length) for _ in range(m)]
    # Iterate through the pattern
    for i in range(m):
        # Iterate through the alphabet
        for j in range(alphabet_length):
            # If the character is an unknown character (wildcard character)
            if pat[i] == unknown_char and bad_char_table[i - 1][j] != -1:
                # Update all indices of rightmost occurences of characters in the pattern
                bad_char_table[i][j] = i
            else:
                # Copy the rightmost occurence of the character from the previous iteration
                bad_char_table[i][j] = bad_char_table[i - 1][j]
            # If the ord value of the character matches the index of the alphabet
            if (ord(pat[i]) - alphabet_start) == j:
                # Store the index of the rightmost occurence of the character 
                bad_char_table[i][j] = i
    return bad_char_table


def boyer_moore(txt, pat):
    """Boyer-Moore pattern matching algorithm.
    
    Performs pattern matching on the text based on a given pattern that
    may contain 0 or more unknown characters (wildcard characters) in the pattern 
    using the edited Boyer-Moore algorithm.

    Args:
        txt (str): The text to search for the pattern.
        pat (str): The pattern to search for in the text.
    
    Returns:
        Prints the starting indices of the occurences of the pattern in the text 
        to an output text file named "output_q2.txt". (1-based indexing)
    """
    global unknown_char, alphabet_start, alphabet_length
    n = len(txt)
    m = len(pat)

    # Preprocessing phase
    bad_char_table = compute_bad_char_table(pat)
    z_suffix = compute_Z_suffix_values(pat)
    good_suffix = compute_good_suffix(z_suffix)
    matched_prefix = compute_matched_prefix(z_suffix, m)

    # Searching phase
    result = []
    start, stop = -1, -1
    # Cumulative shift value of the pattern with respect to the text
    shift = 0
    while shift <= n - m:
        # Set k to the length of the pattern
        k = m - 1
        # Perform explicit comparison of the pattern and the text from right to left
        while k >= 0 and (pat[k] == txt[k + shift] or pat[k] == unknown_char):
            # Skip over the region that is guaranteed to match in the next iteration
            if stop != -1 and (k - 1) == stop:
                k = start - 1
                # Reset the start and stop values
                start, stop = -1, -1
            else:
                # Decrement k as long as the characters match
                k -= 1

        # An occurence of the pattern is found in the text if k < 0
        if k < 0:
            # Append the starting index of the pattern in the text to the results list
            result.append(shift)
            # Shift the pattern to the right by the length of the pattern - matched_prefix[1]
            shift += m - matched_prefix[1]
        else:
            # Determine the shift value based on the good suffix rule
            # If there is no occurence of the good suffix elsewhere in the pattern
            if good_suffix[k + 1] == -1:
                gs_shift = m - matched_prefix[k + 1]
                # Update start and stop values
                start = 1
                stop = matched_prefix[k + 1]
            else:
                # There exist an occurence of the good suffix elsewhere in the pattern
                gs_shift = m - good_suffix[k + 1] - 1
                # Update the start and stop values
                start = good_suffix[k + 1] - m + k + 1
                stop = good_suffix[k + 1]

            # Determine the shift value based on the extended bad character rule
            ebc_shift = k - bad_char_table[k][ord(txt[k + shift]) - alphabet_start]

            # Update the actual shift value based on the maximum of the two shift values
            shift += max(ebc_shift, gs_shift)

    return result


if __name__ == "__main__":
    # Test case 1
    # txt = "ddedadudadededududumadgaergadbvadgaegsdg"
    # pat = "de!!du!"

    # Retrieve the file paths from the command line arguments
    _, txt_file, pat_file = sys.argv

    # Read and store the contents of the text and pattern files
    txt_file = read_file_to_string(txt_file)
    pat_file = read_file_to_string(pat_file)

    # Perform pattern matching using the Boyer-Moore algorithm
    result = boyer_moore(txt_file, pat_file)

    # Print the starting indices of the pattern in the text to an output file
    write_result_to_file("output_q2.txt", result)
