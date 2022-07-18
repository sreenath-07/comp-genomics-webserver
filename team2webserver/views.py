from django.http import HttpResponse, FileResponse, Http404
from django.shortcuts import render, redirect
from django.urls import reverse
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.core.mail import EmailMessage
from django.core.validators import validate_email
from django.template.loader import render_to_string
import os 
import pathlib
import sys
import uuid

import shutil
import os
from subprocess import PIPE, run
from sys import exit


source = os.getcwd() + '/media/'
target = os.getcwd() + '/tmp/'
# output_path = 'output_from_pipeline.txt' # for testing

name='asd'
emailid='asd@gmail.com'

# output_path = ""
global output_path


#output_path = "/home/Team5/Team2-WebServer/Backend/output/sample_userid/output.tar.gz"

def move_to_new_folder(source=source, target=target):
    all_files = os.listdir(source)
    for file in all_files:
        os.rename(source+file, target+file)

def clean_folder(dir=target):
    for f in os.listdir(dir):
        os.remove(os.path.join(dir, f))

# Home Page
def index(request):
    global output_path
    # File upload
    if request.method == 'POST' and 'file_form' in request.POST:
        myfile = request.FILES['myfile']
        fs = FileSystemStorage()
        filename = fs.save(myfile.name, myfile)
        uploaded_file_url = fs.url(filename)
        
        clean_folder()
        move_to_new_folder()

        return render(request, 'index.html', {
            'uploaded_file_url': uploaded_file_url,
        })

    # Email is given
    if request.method == 'POST' and 'username_form' in request.POST:
        email_address = request.POST['emailid']
        invalid = 'yes'

        # Check format validation
        # If valid
        try: 
            validate_email(email_address)
            # redirect('view-results', pk=str(email_address))

            # save user input in text (under same level as manage.py)
            # another way: os.mkdir
            with open(os.getcwd() + '/user_email.txt', 'a') as f:
                f.write(f'{email_address}\n')

            #################################
            # place the pipeline script here 
            # fake_test is just for testing, will write 20*2 in "os.getcwd() + '/output_from_pipeline.txt'"
            #fake_test.multiply(20) 

            name='Team2'
            emailid = email_address
            uniqueid=uuid.uuid1()
        
            current_dir=os.getcwd()
            print(current_dir)
        
            inputpath="/projects/groupb/Team2-WebServer/tmp"
            actuallocation="/projects/groupb/Team2-WebServer/sample_inputs/"

            #create a new id
            new_dictionary={uniqueid.hex: [name,emailid]}
            #print(new_dictionary)
            with open("/projects/groupb/Team2-WebServer/users.log","a+") as log:
              log.write(str(new_dictionary)+"\n")
            current_dir=os.getcwd()
            
            for key in new_dictionary:
              key=key
            user_folder=actuallocation+key
            shutil.copytree(str(inputpath),str(user_folder),copy_function=shutil.move)
            #unzip the files
            
            os.chdir(user_folder)
            path1=os.getcwd()
            
            new_folder=user_folder+'/input'
            
            
            #### cretae input folder and passing it
            os.system("mkdir input")
            os.system('cp -r *.tar.gz input/.')
            
            
            print(new_folder)
            os.chdir(new_folder)
            
            os.system("tar -xvf *.tar.gz")
            os.system("rm -r *.tar.gz")
            
            current_dir=os.getcwd()
            print(current_dir)


            main_script = run([sys.executable, '/projects/groupb/Team2-WebServer/Backend/main.py', '-input_path',user_folder], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            print(main_script.returncode, main_script.stdout, main_script.stderr)
            download_link=path1+'/final_results/final_results.tar.gz'
            output_path=download_link ### download path
            
            
            print(download_link,'\n') #####path of the output folder
            print("Output path of ",output_path)

            ################################

        # If not valid
        except: 
            return render(request, 'index.html', {
                'validation': invalid,
                })

        # send an email with attachment
        mail = EmailMessage(
            'No reply - Analysis result from predictive server.', # subject
            'Greetings from Team2\'s Ecoli Predictive Webserver. \n Your result is attached. \n To unzip your file please use command: tar -C <path_to_untar> -zxvf results.tar.gz \n Inside the tar.gz file you will find fasta files and results from each part of the pipeline along with visualizations. \n Thank you!', # message
            'team2webserver@gmail.com', # from email (Gmail account)
            [email_address], # to email
            )

        # Attach one zip file in specific result folder
        mail.attach_file(output_path)
        mail.send()

        return redirect('view-results', pk=str(email_address))

    return render(request, 'index.html')


# Download Files
def download(request):
    print(output_path) ### download path
    file_server = pathlib.Path(output_path)

    if not file_server.exists():
        return render(request, 'view_results.html', {
                'message': 'file not found.',
            })
    else:
        file_to_download = open(str(file_server), 'rb')
        response = FileResponse(file_to_download, content_type='application/force-download') 
        response['Content-Disposition'] = 'cdattachment; filename="results.tar.gz' # Align with input file type
        return response


# Result Page
def view_results(request, pk):
    with open('/projects/groupb/Team2-WebServer/user_email.txt') as f:
        if pk in f.read():
            return render(request, 'view_results.html', {'email':pk})
        else:
            raise Http404("Page does not exist")



