class PyDAIRUtils:
    """Utilis for PyDAIR application.
    
    Contains some convinient methods that are used in PyDAIR application.
    """
    
    def dot_to_none(self, data):
        """Change `.` to `None` object.
        
        Args:
            data (str, list):
                a string or list class object, that may contain `.`.
        
        Returns:
            The input `data` object is returned. The elements of `.` in `data`
            has been changed to `None` object.
        
        Change all `.` to `None` object in the given a string or list class object.
        The `.` is recoreded in PYDAIR format flat file, therefore, the `.` is 
        converted to `None` object during parsing flat file.
        """
        
        if isinstance(data, list):
            for i in range(len(data)):
                if data[i] == '.':
                    data[i] = None
        else:
            if data == '.':
                data = None
        return data
    
    
    
    
    
    def none_to_dot(self, data, convert_to_str = False, empty_to_dot = True):
        """Changes `None` object to `` (null character) or `.`.
        
        Args:
            data (str, list):
                a string or a list class object that may contain `None` object.
            convert_to_str (bool):
                If `True`, then convert object to string.
            empty_to_dot (bool):
                If `True`, then convert  `` to `.`. Default is `True`.
        
        Returns:
            The input `data` object is returned. The elements of `None`
            has been changed to `.` object.
        
        Change all `None` object to `.` in the given a string or list class object.
        The `data` may contains the object that is not string or None, at that time,
        specify `convert_to_str = True` to forcibly convert object into string.
        """
        
        if isinstance(data, list):
            for i in range(len(data)):
                if data[i] is None:
                    data[i] = '.'
                if convert_to_str:
                    data[i] = str(data[i])
                if empty_to_dot and data[i] == '':
                    data[i] = '.'
        else:
            if data is None:
                data = '.'
            if convert_to_str:
                data = str(data)
            if empty_to_dot and data == '':
                data = '.'
        return data



