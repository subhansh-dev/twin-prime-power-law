"""
Segmented sieve to count twin primes up to 10^10.
Uses block processing to handle the large range without excessive memory.
"""
import math

def simple_sieve(limit):
    """Standard sieve of Eratosthenes up to limit."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if is_prime[i]:
            for j in range(i*i, limit + 1, i):
                is_prime[j] = False
    return is_prime

def count_twins_in_block(block_start, block_size, is_prime_small, primes_small):
    """Count twin primes in [block_start, block_start + block_size)."""
    block = [True] * block_size
    for p in primes_small:
        if p * p > block_start + block_size:
            break
        # Find the first multiple of p in the block
        start = max(p * p, block_start)
        start = start + (p - start % p) % p
        if start == block_start + block_size:
            continue
        if start % p == 0 and start != p:
            block[start - block_start] = False
        for j in range(start + p, block_start + block_size, p):
            block[j - block_start] = False
    
    count = 0
    for i in range(block_size - 2):
        if block[i] and block[i + 2]:
            count += 1
    return count

def count_twin_primes(limit):
    """Count twin primes up to limit using segmented sieve."""
    # Step 1: Simple sieve up to sqrt(limit)
    sqrt_limit = int(math.isqrt(limit)) + 1
    is_prime_small = simple_sieve(sqrt_limit)
    primes_small = [i for i in range(2, sqrt_limit + 1) if is_prime_small[i]]
    
    # Step 2: Segmented sieve in blocks
    block_size = 10_000_000  # 10M numbers per block
    twin_count = 0
    
    # Handle the small range first (up to sqrt_limit)
    for i in range(sqrt_limit - 1):
        if is_prime_small[i] and is_prime_small[i + 2]:
            twin_count += 1
    
    # Segmented blocks
    for block_start in range(sqrt_limit, limit, block_size):
        actual_block_size = min(block_size, limit - block_start)
        count = count_twins_in_block(block_start, actual_block_size, is_prime_small, primes_small)
        twin_count += count
        
        if block_start % (100 * block_size) == 0:
            print(f"  Processed up to {block_start + actual_block_size:,} / {limit:,}")
    
    return twin_count

if __name__ == "__main__":
    limit = 10**10
    print(f"Counting twin primes up to {limit:,}...")
    result = count_twin_primes(limit)
    print(f"pi_2({limit:,}) = {result:,}")
    print(f"\nNicely's verified value: 27,412,679")
    print(f"Match: {result == 27412679}")
