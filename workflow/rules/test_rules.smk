rule test_notebook:
    input:
        'test_{cline}.txt'
    output:
        'test_{cline}.out'
    log:
        notebook="test_notebook_{cline}.ipynb"
    resources:
        tmpdir="./test_scratch"
    notebook:
        "../scripts/sgls/test_notebook.py.ipynb"
