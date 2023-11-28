# Configure Visual Code Studio
Goal: edit python, bash, ... using VSC on mac connected to gi private server (mustard, ...)

```
Andrew E. Davidson
aedavids@ucsc.edu
9/13/23
```

**0) Prerequisits** 
- you can ssh to remote server
- you have cloned the repo
- you have installed the required conda enviroment

**1 ) Install Required Extensions:** 

- Python: Provides Python language support, debugging, and other essential features.

- Visual Studio Code Remote - SSH: Allows you to connect to a remote server via SSH and run VSCode as if it were local.

- Visual Studio Code Remote - SSH: Editing: This is a companion extension for Remote SSH.
To install these extensions:

Open VSCode.
Go to the Extensions view by clicking on the square icon on the sidebar or press Ctrl+Shift+X.
Search for the above extensions and install them.

** 2) Set up Remote SSH:**  
ref: [ see 'Connecting using SSH' https://code.visualstudio.com/docs/remote/ssh-tutorial](https://code.visualstudio.com/docs/remote/ssh-tutorial)

- use command 'Connect To Host'
- on the welcome tab under start click on Open

** 3) Specify Conda Enviroment In VSCode
- Open the terminal in VSCode (Ctrl+ or View > Terminal).

- activate conda env
  ```
  aedavids@mustard $ conda activate extraCellularRNA
  ```
1. With the environment activated, in the bottom-left of VSCode, you should see the Python version. Click on it.

2. This will prompt a list of available interpreters. Choose the one that corresponds to your Conda environment.

3. If it's not listed, you might need to set the path manually:
   * Open the command palette with Ctrl+Shift+P.
   * Search for "Python: Select Interpreter".
   * Click on "Enter interpreter path" and provide the path to the Python executable in your Conda environment. This is typically located in the bin (or Scripts on Windows) directory of the environment folder.

**3) Configure PYTHON PATH**

There are several ways you can set enviroment variable. You can use a extraCellularRNA/.vscode/lauch.json, extraCellularRNA/.vscode/setting.json or extraCellularRNA/.env

launch.json did not work for me. I was unable to run unit test. imports of local python failed. I think the issue is we need to define PYTHONPATH in vscode process before the lanch.json script is run? I think you might need launch.json to configure debugging?  extraCellularRNA/.env did not work. Here is a ref for setting up the pythonpath using [launch.setting](https://k0nze.dev/posts/python-relative-imports-vscode/)



extraCellularRNA/.vscode/setting.json works. see [How can I globally set the PATH environment variable in VS Code?](https://stackoverflow.com/a/45637716).

Here is an example of my settings.json
```
aedavids@mustard $ cat extraCellularRNA/.vscode/settings.json 
{
    "terminal.integrated.env.linux": {
        "PYTHONPATH":"/private/home/aedavids/extraCellularRNA/deconvolutionAnalysis/python"
    },
    "terminal.integrated.env.osx": {
        "PYTHONPATH":"/Users/andrewdavidson/googleUCSC/kimLab/extraCellularRNA/deconvolutionAnalysis/python"
    }
}
```
