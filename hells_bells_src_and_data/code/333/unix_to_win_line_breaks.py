import os

directory = 'successful_algorithms'

for filename in os.listdir(directory):
    if filename.endswith('.txt'):  # Change this to your file type
        file_path = os.path.join(directory, filename)
        with open(file_path, 'rb') as file:
            content = bytearray(file.read())
        content = content.replace(b'\x0D\x0A', b'\x0A')
        with open(file_path, 'wb') as file:
            file.write(content)