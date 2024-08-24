unknown_char = "!"
alphabet_start = 37 

def compute_Z_suffix_values(pat):
    m = len(pat)
    global unknown_char
    # Initialise an array of Z values of size m
    Z_values = [0] * m
    # Reverse the pattern
    reverse_pat = pat[::-1]
    #print(reverse_pat)
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
            while i < m and (reverse_pat[i] == reverse_pat[i-k] or reverse_pat[i] == unknown_char or reverse_pat[i - k] == unknown_char):

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
                while i < m and (reverse_pat[i] == reverse_pat[i-k] or reverse_pat[i] == unknown_char or reverse_pat[i - k] == unknown_char):
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
    n = len(Z_suffix)
    print(n)
    good_suffix = [-1] * (n+1)

    for p in range(n):
        print(p)
        j = n - Z_suffix[p]
        good_suffix[j] = p

    return good_suffix

def compute_matched_prefix(Z_suffix, n):

    matched_prefix = [0] * (n + 1)

    j = 0
    for k in range(n-1, 0, -1):
        if k + Z_suffix[j] == n:
            matched_prefix[k] = Z_suffix[j]
            j+=1
        else:
            matched_prefix[k] = matched_prefix[k+1]
            j += 1
    matched_prefix[0] = n
    return matched_prefix

def boyer_moore(txt, pat):
    n = len(txt)
    m = len(pat)
    global unknown_char
    # Preprocessing phase
    z_suffix = compute_Z_suffix_values(pat)
    good_suffix = compute_good_suffix(z_suffix)
    matched_prefix = compute_matched_prefix(z_suffix, m)

    print(z_suffix)
    print(good_suffix)
    print(matched_prefix) 

    # Searching phase
    result = []
    start, stop = -1, -1
    # Cumulative shift value of the pattern with respect to the text
    shift = 0
    while shift <= n - m:
        # Set i to the length of the pattern
        k = m - 1
        #print(shift)
        # Perform explicit comparison of the pattern and the text from right to left
        while k >= 0 and (pat[k] == txt[k + shift] or pat[k] == unknown_char):
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
            if good_suffix[k+1] == -1:
                print("here")
                gs_shift = m - matched_prefix[k+1]
                # Update start and stop values
                start = 1
                stop = matched_prefix[k+1]
            else:
                gs_shift = m - good_suffix[k+1]
                # Update the start and stop values
                start = good_suffix[k+1] - m + k + 1
                stop = good_suffix[k+1]
            # Update the shift value
            shift += gs_shift
    
    return result

if __name__ == "__main__":
    # Test case 1
    txt = "ddedadudadededududumadgaergadbvadgaegsdg"
    pat = "de!!dud"

    print(boyer_moore(txt, pat))

    # Z_suffix = compute_Z_suffix_values(pat)
    # print(Z_suffix)
    # good_suffix = compute_good_suffix(Z_suffix)
    # print(good_suffix)
    # matched_prefix = compute_matched_prefix(Z_suffix, len(pat))
    # print(matched_prefix)