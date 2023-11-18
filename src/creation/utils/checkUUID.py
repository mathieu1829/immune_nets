import uuid

def checkUUID(id):
    try:
        other_id = uuid.UUID(id, version = 4)
    except:
        return False
    return id == other_id.hex
