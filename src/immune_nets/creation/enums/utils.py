def makeEnumDict(enum):
    dict = enum.__dict__
    return {key:dict[key] for key in dict if not key.startswith('__') and not callable(getattr(enum, key))}

def getEnumValueList(enum):
    dict = enum.__dict__
    return [dict[key] for key in dict if not key.startswith('__') and not callable(getattr(enum, key))]




