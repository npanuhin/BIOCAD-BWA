from functools import cmp_to_key


def comp(s_int, x, y):
    n = len(s_int)
    for i in range(n):
        res = s_int[(x + i) % n] - s_int[(y + i) % n]
        if res < 0:
            return -1
        if res > 0:
            return 1
    return 0


def bwt(s):
    s = s + "$"
    n = len(s)

    s_int = [ord(s[i]) for i in range(n)]

    starts = sorted(list(range(n)), key=cmp_to_key(lambda x, y: comp(s_int, x, y)))

    return "".join([s[(starts[i] + n - 1) % n] for i in range(n)])


def ibwt(s):

    def gen(k):
        permutation = sorted((t, i) for i, t in enumerate(s))
        for _ in range(len(s) - 1):
            t, k = permutation[k]
            yield t

    return ''.join(gen(s.find('$')))


def compress(s):
    res = str(s[0])
    count = 1
    for i in range(1, len(s)):
        if s[i - 1] != s[i]:
            if count > 1:
                res += str(count)
            count = 0
            res += s[i]
        count += 1

    if count > 1:
        res += str(count)

    return res


if __name__ == "__main__":

    # text = input("Input: ")
    with open("src/BWT_test.txt", 'r', encoding="utf-8") as file:
        text = file.read().strip().replace("\n", ' ')[:100000]

    print("Encrypting...")

    encrypted_text = bwt(text)

    print("Encrypted")
    print("Compressing...")

    compressed_text = compress(encrypted_text)

    print("Compressed")

    decrypted_text = ibwt(encrypted_text)

    print("\nText:")
    print(text)
    print("\nEncrypted_text:")
    print(encrypted_text)
    print("\nCompressed_text:")
    print(compressed_text)
    print("\nDecrypted_text:")
    print(decrypted_text)

    # print(len(text), "->", len(encrypted_text), "->", len(compressed_text))

    compress_raw_text = compress(text)
    # print()
    # print(compress_raw_text)
    print(len(text), "->", len(encrypted_text), "->", len(compressed_text), "/", len(compress_raw_text))
