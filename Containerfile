FROM fedora:39
MAINTAINER Gerardo Zapata gerardo.zapata@mcgill.ca 

RUN mkdir /app 
WORKDIR /app

# Download python 
RUN dnf install -y python3.10
RUN dnf install -y python3-pip
COPY requirements.txt .
RUN python3 -m pip install --no-cache-dir --upgrade -r requirements.txt

# Copy the app
COPY PARPAL.py .

RUN chmod 755 PARPAL.py

## from project_tracking
EXPOSE 8000

# Run app on port 8080
ENTRYPOINT ["uvicorn", "PARPAL:app", "--host", "0.0.0.0", "--port", "8000"]
