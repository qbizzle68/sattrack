import json


class ComplexEncoder(json.JSONEncoder):
    """Use this class as the 'cls' argument in the 'json.dumps' method to enable serialization of supported sattrack types."""
    def default(self, obj):
        if hasattr(obj, 'toJson'):
            return obj.toJson()
        else:
            return json.JSONEncoder.default(self, obj)