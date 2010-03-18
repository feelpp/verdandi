#define VERDANDI_LOG_IS_ACTIVE false
#define VERDANDI_LOG_FILENAME "verdandi-%{D}.log"

#include "Verdandi.hxx"

using namespace Verdandi;

#include "Logger.cxx"


class ClassTest
{
public:
    string GetName() const
    {
        return "ClassTest";
    }
    void MemberFunction()
    {
        Logger::Log(*this, "ok");
    }
};


int main(int argc, char** argv)
{
    TRY;

    Logger::Log<5>("TEST 1", "ok");

    Logger::Activate();

    Logger::Log<5>("TEST 2", "ok");

    Logger::SetOption(Logger::stdout_ | Logger::file_, true);

    Logger::Log<5>("TEST 3", "ok");

    Logger::Log<-5>("TEST 4", "ok");

    Logger::Command("hline", "-", Logger::file_);

    Logger::InitializeOptions();

    ClassTest test;
    test.MemberFunction();

    END;

    return 0;
}
