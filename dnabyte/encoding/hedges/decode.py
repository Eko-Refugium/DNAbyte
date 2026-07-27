import traceback


from dnabyte.encoding.hedges.hedges_python import decode_hedges


def decode(data, params, logger=None):
    """
    Decodes HEDGES DNA strands back into the original binary string.
    """

    try:
        # data.data is expected to be a list of DNA strings
        decoded_bytes = decode_hedges(data.data)

        # Convert bytes back to a binary string
        decoded_binary = ''.join(f'{b:08b}' for b in decoded_bytes)

        print("====================================================================")
        print(decoded_bytes)
        print(decoded_binary)
        print("====================================================================")

        info = {
            "number_of_strands": len(data.data),
            "strand_length": len(data.data[0]) if len(data.data) else 0,
            "decoded_bytes": len(decoded_bytes),
        }

        return decoded_binary, True, info

    except Exception as e:
        if logger:
            logger.error(f"Error during HEDGES decoding: {e}")
            logger.error(traceback.format_exc())

        return None, False, {}