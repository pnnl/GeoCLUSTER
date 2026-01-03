import os
from app import app

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8060))
    app.run(
        port=port,
        debug=True,
        ssl_context="adhoc",
    )

