FROM python:3.10.8

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . ./
RUN chmod +x /app/GC_percentage_masked_unmasked_in_genome.py

# Install any needed Python packages
RUN pip install argparse biopython matplotlib pandas numpy

# Make port 80 available to the world outside this container
EXPOSE 8080

# Run your Python script when the container launches
ENTRYPOINT ["python", "GC_percentage_masked_unmasked_in_genome.py"]
