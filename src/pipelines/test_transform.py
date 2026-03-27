from src.services.transformation_service import TransformationService

if __name__ == "__main__":
    service = TransformationService()
    df = service.transform()
    print(df.head())