'''Main entry point for getphylo.'''

import logging
from datetime import datetime
from getphylo import main
from getphylo.utils.errors import GetphyloError

def entrypoint():
    '''Entry point for getphylo'''
    main.initialize_logging()
    try:
        start_time = datetime.now()
        main.main()
        end_time = datetime.now()
        run_time = end_time - start_time
        logging.info('Thank you for using getphylo. The analysis took %s', run_time)
    except GetphyloError as error:
        logging.error(error)
        exit(1)
