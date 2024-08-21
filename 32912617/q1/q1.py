def compute_Z_suffix_values(pat):
    
    n = len(pat)
    Z_values = [-1] * n

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

    matched_prefix = [-1] * (n + 1)
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


if __name__ == "__main__":
    txt = "acababacabacaba"
    pat = "acababacaba"
    print(len(pat))
    z_suffix = compute_Z_suffix_values(pat)
    print(z_suffix)
    print(compute_good_suffix(z_suffix))
    print(compute_matched_prefix(z_suffix, len(pat)))