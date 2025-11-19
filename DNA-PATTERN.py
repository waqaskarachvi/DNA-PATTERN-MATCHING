import streamlit as st
import pandas as pd
from collections import deque, defaultdict

# -----------------------------------------------------------------------------------
# SEARCH ALGORITHMS (SINGLE PATTERN)
# -----------------------------------------------------------------------------------

def naive_search(text, pattern):
    occurrences = []
    n = len(text)
    m = len(pattern)
    for i in range(n - m + 1):
        if text[i:i + m] == pattern:
            occurrences.append(i)
    return occurrences


def compute_lps(pattern):
    lps = [0] * len(pattern)
    j = 0
    i = 1
    while i < len(pattern):
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j
            i += 1
        else:
            if j != 0:
                j = lps[j - 1]
            else:
                lps[i] = 0
                i += 1
    return lps


def kmp_search(text, pattern):
    lps = compute_lps(pattern)
    occurrences = []
    i = j = 0
    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1
        if j == len(pattern):
            occurrences.append(i - j)
            j = lps[j - 1]
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return occurrences


def bad_character_table(pattern):
    table = {}
    length = len(pattern)
    for i in range(length):
        table[pattern[i]] = i
    return table


def boyer_moore_search(text, pattern):
    occurrences = []
    bad_char = bad_character_table(pattern)
    n = len(text)
    m = len(pattern)
    s = 0
    while s <= n - m:
        j = m - 1
        while j >= 0 and pattern[j] == text[s + j]:
            j -= 1
        if j < 0:
            occurrences.append(s)
            s += (m - bad_char.get(text[s + m - 1], -1)) if s + m < n else 1
        else:
            s += max(1, j - bad_char.get(text[s + j], -1))
    return occurrences


def rabin_karp_search(text, pattern, d=256, q=101):
    n = len(text)
    m = len(pattern)
    occurrences = []
    h = pow(d, m - 1) % q
    p = 0
    t = 0
    for i in range(m):
        p = (d * p + ord(pattern[i])) % q
        t = (d * t + ord(text[i])) % q
    for s in range(n - m + 1):
        if p == t:
            if text[s:s + m] == pattern:
                occurrences.append(s)
        if s < n - m:
            t = (d * (t - ord(text[s]) * h) + ord(text[s + 1 + m - 1])) % q
            if t < 0:
                t += q
    return occurrences


# -----------------------------------------------------------------------------------
# AHO–CORASICK (MULTI-PATTERN)
# -----------------------------------------------------------------------------------

class AhoNode:
    def __init__(self):
        self.children = {}
        self.fail = None
        self.output = []


class AhoCorasick:
    def __init__(self, patterns):
        self.root = AhoNode()
        self.build_trie(patterns)
        self.build_failure_links()

    def build_trie(self, patterns):
        for index, pat in enumerate(patterns):
            node = self.root
            for char in pat:
                if char not in node.children:
                    node.children[char] = AhoNode()
                node = node.children[char]
            node.output.append((index, pat))

    def build_failure_links(self):
        queue = deque()
        for child in self.root.children.values():
            child.fail = self.root
            queue.append(child)

        while queue:
            current = queue.popleft()
            for char, next_node in current.children.items():
                fail_node = current.fail
                while fail_node and char not in fail_node.children:
                    fail_node = fail_node.fail
                next_node.fail = fail_node.children[char] if fail_node and char in fail_node.children else self.root
                next_node.output += next_node.fail.output
                queue.append(next_node)

    def search(self, text):
        node = self.root
        results = defaultdict(list)
