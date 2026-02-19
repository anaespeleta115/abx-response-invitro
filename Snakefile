rule hello:
    output:
        "test.txt"
    shell:
        "echo hello > test.txt"
