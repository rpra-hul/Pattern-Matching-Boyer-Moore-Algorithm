def compute_Z_suffix_values(pat):
    
    n = len(pat)
    Z_values = [0] * n

    reverse_pat = pat[::-1]
    l, r, i = 0, 0, 0
    for k in range(1, n):
        # Base Case: k = 1 or k > r, Perform explicit character comparison
        if k > r:
            i = k
            count = 0
            while i < n and reverse_pat[i] == reverse_pat[i-k]:
                i += 1
                count += 1
            Z_values[k] = count
            l = k
            r = k - 1 + count
        # Case 2: k <= r, Perform implicit character comparison
        elif k <= r:
            rem = r - k + 1
            k_prime = k - l
            if Z_values[k_prime] < rem:
                Z_values[k] = Z_values[k_prime]
            elif Z_values[k_prime] == rem:
                i = r
                count = 0
                while i < n and reverse_pat[i] == reverse_pat[i-k]:
                    i += 1
                    count += 1
                Z_values[k] = count
                l = k
                r = k - 1 + count
            else:
                # Z[k'] > rem
                Z_values[k] = rem

    Z_suffix = Z_values[::-1]
    return Z_suffix

def compute_good_suffix(Z_suffix):
    n = len(Z_suffix)
    good_suffix = [-1] * (n+1)

    for p in range(n - 1):
        j = n - Z_suffix[p]
        if good_suffix[j] == -1:
            good_suffix[j] = [p]    
        else:
            good_suffix[j].append(p)

    return good_suffix

def compute_matched_prefix(Z_suffix, n):

    matched_prefix = [0] * (n + 1)
    #print(n, len(Z_suffix), len(matched_prefix))

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

    # Preprocessing phase
    z_suffix = compute_Z_suffix_values(pat)
    good_suffix = compute_good_suffix(z_suffix)
    matched_prefix = compute_matched_prefix(z_suffix, m)

    # Searching phase
    result = []
    start, stop = -1, -1
    # Cumulative shift value of the pattern with respect to the text
    shift = 0
    while shift < n - m:
        # Set i to the length of the pattern
        k = m - 1
        
        # Perform explicit comparison of the pattern and the text from right to left
        while k <= 0 and pat[k] == txt[k + shift]:
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
                gs_shift = m - matched_prefix[k+1]
            else:
                # If there is an occurence of the suffix in the pattern
                bad_char = txt[k + shift]
                for i in range(len(good_suffix[k+1])):
                    if good_suffix[k+1][i] >= k:
                        gs_shift = m - good_suffix[k+1]
                        break
                    else:
                        gs_value = good_suffix[k+1][i]
                        z_suffix_value = z_suffix[gs_value]
                        if gs_value - z_suffix_value >= 0 and pat[gs_value - z_suffix_value] == bad_char:
                            gs_shift = m - gs_value
                            break
    
    return result
    

if __name__ == "__main__":
    txt = "acababacabacaba"
    pat = "acababacaba"
    print(len(pat))
    z_suffix = compute_Z_suffix_values(pat)
    print(z_suffix)
    print(compute_good_suffix(z_suffix))
    print(compute_matched_prefix(z_suffix, len(pat)))