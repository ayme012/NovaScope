# NovaScope: standard-library Python only, so no pip install step.
FROM python:3.12-slim
WORKDIR /app
COPY backend ./backend
COPY frontend ./frontend
ENV HOST=0.0.0.0 PORT=8000 OPEN_BROWSER=0
EXPOSE 8000
CMD ["python", "backend/app.py", "--no-browser"]
