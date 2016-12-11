class PyDAIRUtils:
    """PyDAIRUtils class.
        
    PyDAIRUtils class contains some convinient methods that are
    used in PyDAIR application.
    """
    
    def dot_to_none(self, data):
        """Change `.` to `None` object.
        
        Args:
            data (str): A string or list class object, that may contain `.`.
        
        Returns:
            The input **data** object is returned. However, the elements of `.`
            in **data** object has been changed to `None` object.
        
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
        """Change `None` object to null character or `.`.
        
        Args:
            data (str): A string or a list class object that may contain `None` object.
            convert_to_str (bool): Default is `False`. If `True`, convert any object
                                   (may be integer or others) to string object.
            empty_to_dot (bool): Default is `True`. If `True`, convert null
                                 character to `.`.
        
        Returns:
            The input **data** object is returned. The elements of `None`
            has been changed to null character (when **empty_to_dot =** `False`)
            or `.` (when **empty_to_dot =** `True`) object.
        
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



