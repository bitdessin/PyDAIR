class PyDAIRUtils:
    def dot_to_none(self, data):
        if isinstance(data, list):
            for i in range(len(data)):
                if data[i] == '.':
                    data[i] = None
        else:
            if data == '.':
                data = None
        return data
    
    def none_to_dot(self, data, convert_to_str = False, empty_to_dot = True):
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



