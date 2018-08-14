How to use Python virtualenv

What is Virtualenv?
A Virtual Environment, put simply, is an isolated working copy of Python which
allows you to work on a specific project without worry of affecting other projects
It enables multiple side-by-side installations of Python, one for each project.
It doesn't actually install separate copies of Python, but it does provide a
clever way to keep different project environments isolated. 

https://virtualenv.pypa.io/en/stable/#


Verify if Virtualenv is installed
Run the following command in your terminal
virtualenv --version



#------- project creation --------------------
#1) virtualenv project creation:
#You can have different python env (projects) as many as you want
#2) make a folder to regroup them. Example
mkdir ~/venv_projects
#example for dadi:
cd ~/venv_projects
virtualenv dadi

#------ python module/software installation for a project --------------
#1) we activate this virtual env
source ~/venv_projects/dadi/bin/activate

#you should see (your prompt) something like:
#(dadi) mylogin@pac-sm-gafl01:~$ 

#2) as usual install your packages:
#default version: pip install numpy
#reinstall version: pip install numpy --force-reinstall
#specific version:
pip install 'numpy==1.9.1'
pip install 'scipy==0.15.1'
pip install 'matplotlib==2.1.1'

#install dadi:
tar -xvzf dadi-1.7.0_modif.tar.gz
cd dadi-1.7.0_modif
python setup.py build -f
python setup.py install -f

#3) when finished deactivate the virtualenv:
deactivate


#--- usage ---------------
#1) we activate the virtual env
source ~/venv_projects/dadi/bin/activate
#2) you can run your script:
python myscript.py ...
#3) when finished deactivate the virtualenv:
deactivate


#voila ...
#jacques
