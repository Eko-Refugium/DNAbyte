import os 

class Data:
    """
    Data is the base class for all data types in the DNAbyte pipeline.

    It can be instantiated 
     - with a list of absolute paths to files or, alternatively
     - with a path to a folder using the alternative constructor from_folder.

    :param file_paths: A list of absolute paths to files.    
    """
    
    # Class constant for maximum file size (1MB in bytes)
    MAX_FILE_SIZE = 1024 * 1024  # 1MB
    
    def __init__(self, file_paths):
        """
        Initialize Data object with file paths.
        
        :param file_paths: List of absolute paths to files
        :raises TypeError: If file_paths is not a list
        :raises ValueError: If file_paths is empty, contains invalid paths, or files exceed size limit
        """
        # Validate input type
        if not isinstance(file_paths, list):
            raise TypeError(f"file_paths must be a list, got {type(file_paths).__name__}")
        
        # Check if list is empty
        if len(file_paths) == 0:
            raise ValueError("file_paths cannot be empty")
        
        # Validate each file path
        for file_path in file_paths:
            if not isinstance(file_path, str):
                raise TypeError(f"All file paths must be strings, got {type(file_path).__name__}")
            if not os.path.isfile(file_path):
                raise ValueError(f"File path {file_path} does not point to an existing file.")
            
            # Check file size
            file_size = os.path.getsize(file_path)
            if file_size > self.MAX_FILE_SIZE:
                raise ValueError(f"File {file_path} exceeds maximum size limit of "
                               f"{self.MAX_FILE_SIZE / (1024*1024):.1f}MB. "
                               f"File size: {file_size / (1024*1024):.2f}MB")
        
        self.file_paths = file_paths
        self.size = self.calculate_total_bytes()

    @classmethod
    def from_folder(cls, folder_path):
        """
        Create Data object from all files in a folder.
        
        :param folder_path: Path to directory containing files
        :return: Data object with all files in the folder
        :raises TypeError: If folder_path is not a string
        :raises ValueError: If folder path doesn't exist or contains no files
        """
        if not isinstance(folder_path, str):
            raise TypeError(f"folder_path must be a string, got {type(folder_path).__name__}")
        
        if not os.path.isdir(folder_path):
            raise ValueError(f"Folder path {folder_path} does not point to an existing directory.")
        
        # Get only files, not directories
        all_items = os.listdir(folder_path)
        file_paths = [os.path.join(folder_path, item) for item in all_items 
                     if os.path.isfile(os.path.join(folder_path, item))]
        
        if not file_paths:
            raise ValueError(f"No files found in directory {folder_path}")
        
        return cls(file_paths)

    def calculate_total_bytes(self):
        """
        Calculate total size of all files in bytes.
        
        :return: Total size in bytes
        """
        total_bytes = 0
        for file_path in self.file_paths:
            total_bytes += os.path.getsize(file_path)
        return total_bytes

    def __str__(self):
        """
        String representation of the Data object.
        
        :return: Formatted string with object information
        """
        output = f"Type: {type(self).__name__}\n"
        output += f"Total size: {self.size} bytes ({self.size / (1024*1024):.2f} MB)\n"
        output += f"Number of files: {len(self.file_paths)}\n"
        output += f"File paths: {self.file_paths}\n"
        return output