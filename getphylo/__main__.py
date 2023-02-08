import logging
from getphylo import main
from getphylo.utils.errors import GetphyloError 

def entrypoint():
    main.initialize_logging()
    try:
        main.main()
    except GetphyloError as error:
        logging.error(error)
        exit(1)